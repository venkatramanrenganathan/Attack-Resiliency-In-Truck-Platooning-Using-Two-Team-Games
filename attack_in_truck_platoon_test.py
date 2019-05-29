# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:16:00 2019

@author: vxr131730
"""
import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

# Global Variables
N   = 10            # Total Number of Trucks except the leader
N_b = 3             # Total Number of Bad Trucks
T_s = 0.5           # Sampling Time
Kc  = 3             # Time 
K   = Kc / T_s      # Scaled Time
T   = K + 5         # Total Simulation Steps

class TruckPlatoon():
    # Truck Node Class Object
    def __init__(self):
        
            
    def planTruckPlatoon(self): 
        
        # Get the continuous time platoon dynamics data
        A_c, B_c, B_k, F_k = self.getPlatoonContinousDynamics()
        
        # Get Best Attacker Position
        attackerPositions, gammaRecord = self.getBestAttackerPosition()
        
        # Get good and bad truck position indices
        posParam                                = [attackerPositions, gammaRecord]
        defender_trucks, attacker_trucks, gamma = self.getTruckIndices(posParam)
        
        # Get the discrete time platoon dynamics data
        discreteParam  = [A_c, B_c, B_k, F_k]
        A, B, B_k, F_k = self.getPlatoonDiscreteDynamics(discreteParam)
        
        # Keep only corresponding defender/attacker entries in control matrices
        manipulateParam = [B_k, F_k, defender_trucks, attacker_trucks]
        B_k, F_k        = self.manipulateControlMatrices(manipulateParam)
        
        # Use Backward Induction to Compute Saddle Point Equilibrium
        # Compute the backward recursion coefficients
        inputParam                     = [A, B, B_k, F_k, gamma]
        phi, theta, lambda_1, lambda_2 = self.computeRecursionCoefficients(inputParam)
        
        # Compute optimal nash equilibrium control law for both good and bad trucks
        # x - array of truck positions across whole time horizon
        # u - array of defender truck controls across whole time horizon
        # v - array of attacker truck controls across whole time horizon        
        inputParam = [A, B_k, F_k, theta, phi, lambda_1, lambda_2]
        x, u, v    = self.computeOptimalFeedbackLaw(inputParam) 
        
        # Plot only corresponding controls 
        u[attacker_trucks,:] = []
        v[defender_trucks,:] = []
        
        # Prepare the output paramters and return it
        outParam = [x, u, v, defender_trucks, attacker_trucks]        
        return outParam
    
        
    def getBestAttackerPosition(self):
        # To be filled by Karthik - Hardcoded as of now
        return np.array([2,3,8])
    
    def getPlatoonContinousDynamics(self):
        # Dynamics of Platooning Trucks (except leader) - Continuous time model
        # System Matrix
        A_c = np.diag(np.ones(2*N-1,1),-1) - np.identity(2*N)         
        # Input Matrix
        B_c = np.block([np.zeros((N,N))],
                       [np.diag(np.ones(N-1,1),-1) - np.identity(N)]) 
        # Defender Input Matrix                   
        B_k = B_c
        # Attacker Input Matrix             
        F_k = B_c
        # Return the continuous time dynamics matrices
        return A_c, B_c, B_k, F_k                  
        
     def manipulateControlMatrices(self, manipulateParam):
         # Unbox the input parameters
         B_k             = manipulateParam[0]
         F_k             = manipulateParam[1]
         defender_trucks = manipulateParam[2]
         attacker_trucks = manipulateParam[3]
         for i in range(defender_trucks):
             F_k[:,defender_trucks[i]] = np.zeros((2*N,1))
         for i in range(attacker_trucks):
             B_k[:,attacker_trucks[i]] = np.zeros((2*N,1)) 
         # return the manipulated B_k and F_k matrices
         return B_k, F_k         

    def getPlatoonDiscreteDynamics(self, discreteParam):
        # Unbox the input paramters
        A_c = discreteParam[0]
        B_c = discreteParam[1]
        B_k = discreteParam[2]
        F_k = discreteParam[3]                  
        # Convert to Discrete time model
        A   = np.identity(2*N) + np.dot(T_s,A_c)
        B   = np.dot(T_s, B_c)
        B_k = np.dot(T_s, B_k)
        F_k = np.dot(T_s, F_k) 
        # Return the discrete dynamics 
        return A, B, B_k, F_k     
    
    def getTruckIndices(self, posParam):
        # Unbox the input paramters
        attackerPositions = posParam[0]
        gammaRecord       = posParam[1]
        truck_array       = numpy.arange(1, N, 1)
        # Out of all best positions, get the top N_b positions
        N_b_truck_indices = attackerPositions[:N_b]
        N_b_gamma_indices = gammaRecord[:N_b]
        # Sort in an ascending order and get the best N_b bad truck positions
        attacker_trucks   = sort(N_b_truck_indices)        
        # We can set gamma even to max_gamma = 100000
        gamma             = np.amax(N_b_gamma_indices) + 1 
        # Get the remaining defending truck positions        
        truck_array[attacker_trucks] = [0]* len(attacker_trucks)
        defender_trucks              = np.nonzero(truck_vector)        
        # Return the attacking and defending truck positions
        return defender_trucks, attacker_trucks, gamma
    
    def computeRecursionCoefficients(self, inputParam):
        # Unbox the input paramters
        A     = inputParam[0]  
        B     = inputParam[1]  
        B_k   = inputParam[2]  
        F_k   = inputParam[3] 
        gamma = inputParam[4] 
        # Cost function parameters        
        rho = 0.01                      # penalty scaling constant
        Qf  = 100 * np.identity(2*N)    # Terminal State Penalty Matrix 
        Q   = np.identity(2*N)          # State Penalty Matrix
        R_u = rho * np.identity(N)      # Defender Control Penalty Matrix
        R_v = rho * np.identity(N)      # Attacker Control Penalty Matrix
        D   = 1                         # Desired Safe inter-vehicular distance
        E   = np.dot(D,np.ones((1,2*N)) # Vector specifying the distance offset
        
        # Data Structures to hold backward recursion coefficients
        QfN = np.zeros((2*N,2*N))
        QN  = np.zeros((2*N,2*N))    
        PN  = np.empty((2*N,2*N,K+1))
        qN  = np.empty((2*N,K+1))
        rN  = np.empty((K+1))    
        
        # Simulation Data Structures
        theta     = np.empty((N,2*N,K))
        phi       = np.empty((N,2*N,K))
        lambda_1  = np.empty((N,K))
        lambda_2  = np.empty((N,K))
        alpha     = np.empty((N,N,K))
        alpha_inv = np.empty((N,N,K))
        beta      = np.empty((N,N,K))
        beta_inv  = np.empty((N,N,K))
        
        # Initialize the backward recursion coefficients
        PN[:,:,K+1] = QfN               # P_T = Q_T
        qN[:,K+1]   = np.zeros((2*N,1)) # q_T = 0
        rN[K+1]     = 0                 # r_T = 0
            
        # Zero Sum Case - Value Function is of the form x^T P_t x + 2 q_t' x + r_t
        # Compute the Recursive Coefficents
        for k in range(K,-1,-1):
            alpha[:,:,k]     = R_u + B_k.T @ PN[:,:,k+1] @ B_k        
            alpha_inv[:,:,k] = LA.inv(alpha[:,:,k])    
            beta[:,:,k]      = -np.dot(math.pow(gamma,2), R_v) + F_k.T @ PN[:,:,k+1] @ F_k
            beta_inv[:,:,k]  = LA.inv(beta[:,:,k])    
            mu_k             = alpha_inv[:,:,k] @ B_k.T @ PN[:,:,k+1] @ F_k @ beta_inv[:,:,k] @ F_k.T
            zeta_k           = beta_inv[:,:,k] @ F_k.T @ PN[:,:,k+1] @ B_k @ alpha_inv[:,:,k] @ B_k.T
            kappa            = np.identity(N) - mu_k @ PN[:,:,k+1] @ B_k
            eta              = np.identity(N) - zeta_k @ PN[:,:,k+1] @ F_k    
            theta[:,:,k]     = LA.inv(kappa) @ ((mu_k @ PN[:,:,k+1] - alpha_inv[:,:,k] @ B_k.T @ PN[:,:,k+1]) @ A)
            phi[:,:,k]       = LA.inv(eta) @ ((zeta_k @ PN[:,:,k+1] - beta_inv[:,:,k] @ F_k.T @ PN[:,:,k+1]) @ A)
            lambda_1[:,k]    = LA.inv(kappa) @ ((mu_k - alpha_inv[:,:,k] * B_k.T) @ qN[:,k+1])                
            lambda_2[:,k]    = LA.inv(eta) @ ((zeta_k - beta_inv[:,:,k] @ F_k.T) @ qN[:,k+1])        
            # Calculation of recursion coefficients    
            PN[:,:,k] = Q + A.T @ PN[:,:,k+1] @ A + theta[:,:,k].T @ alpha[:,:,k] @ theta[:,:,k] + 
                        phi[:,:,k].T @ beta[:,:,k] @ phi[:,:,k] + 2 @ theta[:,:,k].T @ B_k.T @ PN[:,:,k+1] @ A +
                        np.dot(2,phi[:,:,k].T) @ F_k.T @ PN[:,:,k+1] @ A + theta[:,:,k].T @ B_k.T @ PN[:,:,k+1] @ F_k @ phi[:,:,k] + 
                        phi[:,:,k].T @ F_k.T @ PN[:,:,k+1] @ B_k @ theta[:,:,k]
        
            qN[:,k]   = (lambda_1[:,k].T @ alpha[:,:,k] @ theta[:,:,k] + lambda_2[:,k].T @ beta[:,:,k] @ phi[:,:,k] + 
                        lambda_1[:,k].T @ B_k.T @ PN[:,:,k+1] @ (A + F_k @ phi[:,:,k]) + 
                        lambda_2[:,k].T @ F_k.T @ PN[:,:,k+1] @ (A + B_k @ theta[:,:,k]) + 
                        qN[:,k+1].T @ (A + B_k @ theta[:,:,k] + F_k @ phi[:,:,k]) - E).T
        
            rN[k]     = rN(k+1) + lambda_1[:,k].T @ alpha[:,:,k] @ lambda_1[:,k] + lambda_2[:,k].T @ beta[:,:,k] @ lambda_2[:,k] + 
                        lambda_1[:,k].T @ B_k.T @ PN[:,:,k+1] @ F_k @ lambda_2[:,k] + 
                        lambda_2[:,k].T @ F_k.T @ PN[:,:,k+1] @ B_k @ lambda_1[:,k] + 
                        np.dot(2*K,D@D) + np.dot(2,qN(:,k+1).T) @ (B_k @ lambda_1[:,k] + F_k @ lambda_2[:,k])
        
        # Return phi, theta, lambda_1, lambda_2
        return phi, theta, lambda_1, lambda_2
        
                        
    def computeOptimalFeedbackLaw(self, inputParam):
        # Defender input, u_k     = theta * x + lambda_1
        # Attacker input, v_k     = phi   * x + lambda_2
        # State Update,   x_{k+1} = A     * x + B_k * u_k + F_k * v_k
        # Unbox all function parameters
        A        = inputParam[0]
        B        = inputParam[1]
        F        = inputParam[2]        
        theta    = inputParam[3]
        phi      = inputParam[4]
        lambda_1 = inputParam[5]
        lambda_1 = inputParam[6]        
        
        # Data Structures to hold values for computing optimal feedback law
        xN = np.empty((2*N,T))  # Truck State Vector 
        uN = np.zeros((N,T))    # Good Truck Controls
        vN = np.zeros((N,T))    # Bad Truck Controls
        # Initialize Trucks to start at origin with negative velocity
        xN(:,1) = np.block([np.zeros(zeros(N,1))], [np.dot(-0.1, np.ones((N,1))])         
        
        # Compute Optimal Feedback Law
        for t in range(T-1):
            if t <= K:
                uN[:,t] = theta[:,:,t] @ xN[:,t] + lambda_1[:,t] 
                vN[:,t] = phi[:,:,t] @ xN[:,t] + lambda_2[:,t]
            xN[:,t+1] = A @ xN[:,t] + B @ uN[:,t] + F @ vN[:,t]
        
        # Return Truck State Vector, Good and Bad truck control sequences
        return xN, uN, vN 
    
    def getTruckTrajectories(self, x):         
        # Obtain the exact truck trajectories
        w = LA.inv(np.block(np.diag(np.ones(N,1),-1) - np.identity(N+1))) @ 
            np.block([np.zeros((N,T))], [x[N+1:,:]])
        trucks_trajectory = -np.ones((2*T,1)) @ numpy.arange(0, N, 1) + np.cumsum(w.T)        
        # return trucks_trajectory
        return trucks_trajectory
    
    def drawTrajectories(self,plotParam, drawFlag):
        plt.clf()   
        # drawFlag == 1 corresponds to plotting truck trajectories
        if drawFlag == 1:
            # Unbox the corresponding plot parameters
            trucksTrajectory = plotParam[0]
            defender_trucks  = plotParam[1]
            attacker_trucks  = plotParam[2]
            if trucksTrajectory is not None:                
                x_values = np.arange(2*T) 
                # Plot the leader truck trajectory
                leader = plt.plot(x_values, trucksTrajectory[:,0], "k")
                # Plot the defender trucks trajectories
                goodTrucks = plt.plot(x_values, trucksTrajectory[:,defender_trucks+1], "b")
                # Plot the attacker trucks trajectories
                badTrucks = plt.plot(x_values, trucksTrajectory[:,attacker_trucks+1], "r")            
                # Add Legends
                plt.legend(handles=[leader, goodTrucks[0], badTrucks[0]])
                plt.xlabel('time (s)')
                plt.ylabel('truck position')
                plt.grid(True)
                plt.show()
        # drawFlag == 2 corresponds to plotting good and bad trucks controls
        elif drawFlag == 2:
            # Unbox the corresponding plot parameters
            goodTrucksCtrl  = plotParam[0]
            badTrucksCtrl   = plotParam[1]
            defender_trucks = plotParam[2]
            attacker_trucks = plotParam[3]
            if goodTrucksCtrl is not None and badTrucksCtrl is not None:
                x_values = np.arange(T)                
                # Plot Good Truck Controls.
                plt.subplot(2,1,1)                
                plt.plot(x_values, goodTrucksCtrl, "b", label='Defending Trucks Control Sequences')                
                plt.legend(loc=1)
                plt.xlabel('time')
                plt.ylabel('Good Trucks Control (acceleration)')
                # Plot Bad Truck Controls.
                plt.subplot(2,1,2)
                plt.plot(x_values, badTrucksCtrl, "r", label='Attacking Trucks Control Sequences')
                plt.legend(loc=1)         
                plt.xlabel('time')
                plt.ylabel('Bad Trucks Control (acceleration)')                
                plt.grid(True)
                plt.show()
                

def main():    
    # Create The class object
    platoon   = TruckPlatoon() 
    # Plan the truck positions along the trajectory
    planParam = platoon.planTruckPlatoon()
    # Unbox the output paramters
    truckPositions  = planParam[0]
    goodTrucksCtrl  = planParam[1]
    badTrucksCtrl   = planParam[2]
    defender_trucks = planParam[3]
    attacker_trucks = planParam[4]    
    # Prepare the Trucks Trajectories    
    trucksTrajectory = platoon.getTruckTrajectories(truckPositions)
    # Plot the Trucks Trajectories 
    plotParam = [trucksTrajectory, defender_trucks, attacker_trucks]
    platoon.drawTrajectories(plotParam, drawFlag = 1)
    # Plot the Good and Bad Trucks Controls
    plotParam = [goodTrucksCtrl, badTrucksCtrl, defender_trucks, attacker_trucks]
    platoon.drawTrajectories(plotParam, drawFlag = 2)
    

if __name__ == '__main__':
    main()