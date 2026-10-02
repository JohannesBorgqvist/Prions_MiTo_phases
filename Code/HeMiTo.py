
# =========================================================
# Script: HeMiTo
# Written by: Johannes Borgqvist
# Date: 2026-10-02
# Description:
# This is the main script for the  continuation paper for the first HeMiTo-paper, and here we validate analytical approximations in the Mi-phase.
# =========================================================
# Import our beloved libraries
from load_libraries import * 
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True
# =========================================================
# =========================================================
# Define dimensionless parameters
# =========================================================
# =========================================================
# Define dimensionless parameters
c1 = 1.75 # Healthy formation rate
c2 = 0.3 # Healthy degradation rate
vareps = 0.05 # Conversion rate and perturbation parameter (small value)
#vareps = 0.10 # Exponential
# =========================================================
# =========================================================
# DEFINE FUNCTIONS
# =========================================================
# =========================================================
# ----------------------------------------------------------
# Function 1: ODE system with f=1
def aB(y, t, p):
    # Unpack and re-name our parameters so that
    # they are understandable
    c1 = p[0]
    c2 = p[1]
    vareps = p[2]
    # Define the system of ODEs
    dydt = [c1-c2*y[0] - vareps*y[0]*y[1], vareps*(y[0]*y[1]-y[1])]
    return dydt
# ----------------------------------------------------------   
# Function 2: ODE system with Hill-type conversion function
def aB_Hill(y, t, p):
    # Unpack and re-name our parameters so that
    # they are understandable
    c1 = p[0]
    c2 = p[1]
    vareps = p[2]
    # Define our non-linear Hill function
    f = 1.90+((y[1]**4)/(1+(y[1]**5)))
    # Define the system of ODEs
    dydt = [c1-c2*y[0] - vareps*f*y[0]*y[1], vareps*(f*y[0]*y[1]-y[1])]
    return dydt
# ----------------------------------------------------------   
# Function 3: Fast He-phase solution
def u_He(t, p):
    # Unpack and re-name our parameters so that
    # they are understandable
    u0 = p[0]
    c1 = p[1]
    c2 = p[2]
    # Define the system of ODEs
    dydt = np.array([(c1/c2)-((c1/c2)-u0)*np.exp(-c2*t_i) for t_i in t])
    return dydt
# ----------------------------------------------------------   
# Function 4: Slow Mi-phase solution
def u_Mi(t, p):
    # Unpack and re-name our parameters so that
    # they are understandable
    u0 = p[0]
    v0 = p[1]
    c1 = p[2]
    c2 = p[3]
    f0 = p[4]
    # Define various terms
    term_1 = v0*f0*(c1/(c2**2))
    # Define the system of ODEs
    dydt = np.array([term_1*(np.exp(-c2*t_i)-1)+((c1/c2)-u0)*t_i*v0*f0*np.exp(-c2*t_i) for t_i in t])
    return dydt
# ----------------------------------------------------------   
# Function 5: Slow Mi-phase solution
def v_Mi(t, p):
    # Unpack and re-name our parameters so that
    # they are understandable
    u0 = p[0]
    v0 = p[1]
    c1 = p[2]
    c2 = p[3]
    f0 = p[4]
    # Define the system of ODEs
    dydt = np.array([v0*((((c1/c2)*f0)-1)*t_i+(f0/c2)*((c1/c2)-u0)*(np.exp(-c2*t_i)-1)) for t_i in t])
    return dydt
# =========================================================
# =========================================================
# MAIN FUNCTION
# =========================================================
def main():
    # ===================================================================
    # PART 1 OF 2: NUMERICAL SOLUTIONS VS TRUNCATED PERTURBATION SERIES
    # ===================================================================
    # Initial conditions
    y0 = [1.00, 0.05]
    # Parameters yeah
    p = [c1, c2, vareps]
    # Define the age or time if you will
    t = np.linspace(0,6,200)
    #-------------------------------------------------------------------
    # Simulate patient with simple conversion term
    prions = odeint(aB, y0, t, args=(p,))
    # Extract each species for patient
    u, v = prions.T
    # Conversion function
    f0 = 1 # Standard heterodimer model
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # UNCOMMENT THE REGION IF YOU WANT TO PLOT THE SOLUTIONS FOR
    # THE NONLINEAR FUNCTION
    #-------------------------------------------------------------------    
    #-------------------------------------------------------------------    
    # Simulate patient with simple conversion term
    prions = odeint(aB_Hill, y0, t, args=(p,))
    # Extract each species for patient
    u, v = prions.T
    # Conversion function f(v0) for the nonlinear Hill conversion function
    f0 = 1.90+((y0[1]**4)/(1+(y0[1]**5)))
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # FIGURE 1: SOLUTIONS VS TRUNCATED PERTURBATION SERIES
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # The colours for the plots are chosen based on the colourmaps available at https://colorbrewer2.org/
    # Make a 1x2 plot
    f_prions, ax_prions = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
    # Plot the simulated data and the underlying model
    ax_prions[0].plot(t,u,color=(2/256,56/256,88/256),label="$u(\\tau)$",linewidth=3.0)
    ax_prions[0].plot(t,u_He(t, [y0[0], c1, c2])+vareps*u_Mi(t, [y0[0], y0[1], c1, c2, f0]),color=(5/256,112/256,176/256),label="$u_{\\mathrm{He}}(\\tau)+\\varepsilon u_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
    # Set a grid and define a legend
    ax_prions[0].grid()
    ax_prions[0].legend(loc='best',prop={"size":40})
    # Set the x-labels and y-labels
    ax_prions[0].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
    ax_prions[0].set_ylabel(ylabel="Particle concentration",fontsize=40)
    ax_prions[0].set_title(label="Healthy species",fontsize=50)    
    ax_prions[0].xaxis.set_tick_params(labelsize=35)
    ax_prions[0].yaxis.set_tick_params(labelsize=35)
    # Plot the simulated data and the underlying model
    ax_prions[1].plot(t,v,color=(0/256,68/256,27/256),label="$v(\\tau)$",linewidth=3.0)
    ax_prions[1].plot(t,y0[1]*((t+1)/(t+1))+vareps*v_Mi(t, [y0[0], y0[1], c1, c2, f0]),color=(102/256,194/256,164/256),label="$v_{\\mathrm{He}}(\\tau)+\\varepsilon v_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
    # Set a grid and define a legend
    ax_prions[1].grid()
    ax_prions[1].legend(loc='best',prop={"size":40})
    # Set the x-labels and y-labels
    ax_prions[1].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
    ax_prions[1].set_ylabel(ylabel="Particle concenctration",fontsize=40)
    ax_prions[1].set_title(label="Toxic species",fontsize=50) 
    ax_prions[1].xaxis.set_tick_params(labelsize=35)
    ax_prions[1].yaxis.set_tick_params(labelsize=35)
    plt.suptitle("Comparison between truncated perturbation series and numerical solutions when $\\varepsilon="+str(vareps)+"$",fontsize=35)
    f_prions.savefig('../Figures/HeMi_dynamics.png')
    # ===================================================================
    # PART 2 OF 2: EXPONENTIAL APPROXIMATION OF TOXIC SPECIES
    # ===================================================================
    # # Parameters yeah
    p_high = [c1, 5.5*c2, vareps]
    p_medium = [c1, 5.0*c2, vareps]
    p_low = [c1, 4.5*c2, vareps]
    # Define the age or time if you will
    t_exp = np.linspace(0,10,200)
    # Simulate patient with simple conversion term
    prions_exp_high = odeint(aB_Hill, y0, t_exp, args=(p_high,))
    prions_exp_medium = odeint(aB_Hill, y0, t_exp, args=(p_medium,))
    prions_exp_low = odeint(aB_Hill, y0, t_exp, args=(p_low,))
    # Extract each species for patient
    u_exp_high, v_exp_high = prions_exp_high.T
    u_exp_medium, v_exp_medium = prions_exp_medium.T
    u_exp_low, v_exp_low = prions_exp_low.T
    # Errors
    err_num_high = np.abs(v_exp_high-y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp))
    err_num_medium = np.abs(v_exp_medium-y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp))
    err_num_low = np.abs(v_exp_low-y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp))
    err_theo_high = vareps*y0[1]*(f0*(c1-5.5*c2*y0[0])-vareps*(y0[0]*f0-1)**2)*((t_exp**2)/2)
    err_theo_medium = vareps*y0[1]*(f0*(c1-5.0*c2*y0[0])-vareps*(y0[0]*f0-1)**2)*((t_exp**2)/2)
    err_theo_low = vareps*y0[1]*(f0*(c1-4.5*c2*y0[0])-vareps*(y0[0]*f0-1)**2)*((t_exp**2)/2)
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # FIGURE 2: SOLUTION OF TOXIC SPECIES VS EXPONENTIAL APPROXIMATION
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # Plot A: Toxic species vs exponential approximation
    #-------------------------------------------------------------------
    #Define the first figure
    f_prions, ax_prions = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
    # Plot the simulated data and the underlying model
    ax_prions[0].plot(t_exp,v_exp_low,color=(35/256,139/256,69/256),label="$v(\\tau)$ low $c_{2}$",linewidth=3.0)
    ax_prions[0].plot(t_exp,v_exp_medium,color=(0/256,109/256,44/256),label="$v(\\tau)$ medium $c_{2}$",linewidth=3.0)
    ax_prions[0].plot(t_exp,v_exp_high,color=(0/256,68/256,27/256),label="$v(\\tau)$ high $c_{2}$",linewidth=3.0)
    ax_prions[0].plot(t_exp,y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp),color=(153/256,216/256,201/256),label="$v_{0}e^{\\lambda\\varepsilon\\tau}$",linewidth=3.0)
    # Set a grid and define a legend
    ax_prions[0].grid()
    ax_prions[0].legend(loc='best',prop={"size":40})
    # Set the x-labels and y-labels
    ax_prions[0].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
    ax_prions[0].set_ylabel(ylabel="Particle concentration",fontsize=40)
    ax_prions[0].set_title(label="Toxic species",fontsize=50)    
    ax_prions[0].xaxis.set_tick_params(labelsize=35)
    ax_prions[0].yaxis.set_tick_params(labelsize=35)
    # Plot the simulated data and the underlying model
    ax_prions[0].plot(t_exp,v_exp_low,color=(35/256,139/256,69/256),label="$v_{\\mathrm{low}}(\\tau)$",linewidth=3.0)
    #-------------------------------------------------------------------
    # Plot B: Error ratio
    #-------------------------------------------------------------------
    ax_prions[1].plot(t_exp,100*np.abs(err_theo_low)/(np.abs(err_theo_low)+np.abs(err_num_low-err_theo_low)),color=(77/256,0/256,75/256),label="$R_{\\mathrm{HSS}}(\\tau)$ low $c_{2}$",linewidth=3.0)    
    ax_prions[1].plot(t_exp,100*np.abs(err_theo_medium)/(np.abs(err_theo_medium)+np.abs(err_num_medium-err_theo_medium)),color=(129/256,15/256,124/256),label="$R_{\\mathrm{HSS}}(\\tau)$ medium $c_{2}$",linewidth=3.0)
    ax_prions[1].plot(t_exp,100*np.abs(err_theo_high)/(np.abs(err_theo_high)+np.abs(err_num_high-err_theo_high)),color=(136/256,65/256,157/256),label="$R_{\\mathrm{HSS}}(\\tau)$ high $c_{2}$",linewidth=3.0)    
    # Set a grid and define a legend
    ax_prions[1].grid()
    ax_prions[1].legend(loc='best',prop={"size":40})
    # Set the x-labels and y-labels
    ax_prions[1].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
    ax_prions[1].set_ylabel(ylabel="$R_{\\mathcal{HSS}}(\\tau)$",fontsize=40)
    ax_prions[1].set_title(label="Error ratio",fontsize=50)    
    ax_prions[1].xaxis.set_tick_params(labelsize=35)
    ax_prions[1].yaxis.set_tick_params(labelsize=35)
    plt.suptitle("Exoponential approximation of the toxic species during the Mi-phase",fontsize=35)
    f_prions.savefig('../Figures/Mi_phase_exponential.png')
# =========================================================
# RUN THE MAIN FUNCTION
# =========================================================
if __name__ == "__main__":
    main()               

