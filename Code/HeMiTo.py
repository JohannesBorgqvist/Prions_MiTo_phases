
# =========================================================
# Script: HeMiTo
# Written by: Johannes Borgqvist
# Date: 2026-09-30
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
    # =========================================================
    # Plot solutions, yeah?
    # =========================================================
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
    # FIGURE 1
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------    
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
    ax_prions[1].plot(t,y0[1]*((t+1)/(t+1))+vareps*v_Mi(t, [y0[0], y0[1], c1, c2, f0]),color=(2/256,56/256,88/256),label="$v_{\\mathrm{He}}(\\tau)+\\varepsilon v_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
    ax_prions[1].plot(t,v,color=(5/256,112/256,176/256),label="$v(\\tau)$",linewidth=3.0)
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
    plt.show()
    f_prions.savefig('../Figures/HeMi_dynamics.png')


# =========================================================
# RUN THE MAIN FUNCTION
# =========================================================
if __name__ == "__main__":
    main()               

# # =========================================================
# # Plot solutions, yeah?
# # =========================================================
# # Initial conditions
# y0 = [1.00, 0.05]
# # Calculate f0=f(v0)
# f0 = 1.9+((y0[1]**4)/(1+(y0[1]**5)))
# print("The initial conversion function value,\t\tf(v0)\t=\t%0.3f\n\n"%(f0))
# #y0 = [0.19, 0.06]
# #y0 = [0.15, 0.10]
# # Parameters yeah
# p = [c1, c2, vareps]
# # Define the age or time if you will
# t = np.linspace(0,6,200)
# #-------------------------------------------------------------------
# # Simulate patient with simple conversion term
# prions_1 = odeint(aB, y0, t, args=(p,))
# # Extract each species for patient
# u_1, v_1 = prions_1.T
# # Simulate patient with simple conversion term
# prions_2 = odeint(aB_Hill, y0, t, args=(p,))
# # Extract each species for patient
# u_2, v_2 = prions_2.T
# #-------------------------------------------------------------------
# # Fast solution
# t2 = np.linspace(0,6,200)
# u_He(t2, [y0[0], c1, c2])
# # Slow solution
# #t3 = np.linspace(0,0.8,50)
# #u_Mi(t3, [y0[0], y0[1], c1, c2, 1.9])
# #---------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------
# # Plot the trajectories in the phase plane
# #---------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------
# # #Define the first figure
# # f_prions, ax_prions = plt.subplots(1, 1, constrained_layout=True, figsize=(20, 8))
# # # #---------------------------------------------------------------------------------
# # # Plot the simulated data and the underlying model
# # ax_prions.plot(t,u_2,color=(2/256,56/256,88/256),label="$u(\\tau)$",linewidth=3.0)
# # ax_prions.plot(t,v_2,color=(5/256,112/256,176/256),label="$v(\\tau)$",linewidth=3.0)
# # ax_prions.plot(np.linspace(0,t[-1],100),np.array([(c1/c2)*(i+2)/(i+2) for i in range(100)]),"--",color=(0,0,0),label="$c_{1}/c_{2}$",linewidth=3)
# # ax_prions.plot(np.array([(1/(y0[0]*1.9-1))*np.log(1/vareps)*(i+2)/(i+2) for i in range(100)]),np.linspace(0,2.5,100),"--",color=(0,0,0),label="$c_{1}/c_{2}$",linewidth=3)
# # # Set a grid and define a legend
# # ax_prions.grid()
# # ax_prions.legend(loc='best',prop={"size":40})
# # # Set the x-labels and y-labels
# # ax_prions.set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
# # ax_prions.set_ylabel(ylabel="Prion conc.",fontsize=40)
# # ax_prions.set_title(label="$f(v)=1.30+\\frac{v^{4}}{1+v^{5}}$",fontsize=50)    
# # ax_prions.xaxis.set_tick_params(labelsize=35)
# # ax_prions.yaxis.set_tick_params(labelsize=35)    
# # f_prions.savefig('./Figures/HeMiTo_dynamics.png')





# #Define the first figure
# f_prions_2, ax_prions_2 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))

# # Plot the simulated data and the underlying model
# ax_prions_2[0].plot(t,u_2,color=(2/256,56/256,88/256),label="$u(\\tau)$",linewidth=3.0)
# ax_prions_2[0].plot(t2,u_He(t2, [y0[0], c1, c2])+vareps*u_Mi(t2, [y0[0], y0[1], c1, c2, f0]),color=(5/256,112/256,176/256),label="$u_{\\mathrm{He}}(\\tau)+\\varepsilon u_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
# #ax_prions_2[0].plot(np.array([(1/(y0[0]*1.9-1))*np.log(1/vareps)*(i+2)/(i+2) for i in range(100)]),np.linspace(0,2.5,100),"--",color=(0,0,0),label="$\\tau=1/\\varepsilon$",linewidth=3)
# # Set a grid and define a legend
# ax_prions_2[0].grid()
# ax_prions_2[0].legend(loc='best',prop={"size":40})
# # Set the x-labels and y-labels
# ax_prions_2[0].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
# ax_prions_2[0].set_ylabel(ylabel="Prion conc.",fontsize=40)
# ax_prions_2[0].set_title(label="Healthy species",fontsize=50)    
# ax_prions_2[0].xaxis.set_tick_params(labelsize=35)
# ax_prions_2[0].yaxis.set_tick_params(labelsize=35)
# # Plot the simulated data and the underlying model
# ax_prions_2[1].plot(t2,y0[1]*((t2+1)/(t2+1))+vareps*v_Mi(t2, [y0[0], y0[1], c1, c2, f0]),color=(2/256,56/256,88/256),label="$v_{\\mathrm{He}}(\\tau)+\\varepsilon v_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
# ax_prions_2[1].plot(t,v_2,color=(5/256,112/256,176/256),label="$v(\\tau)$",linewidth=3.0)
# #ax_prions_2[0][1].plot(t3,u_Mi(t3, [y0[0], y0[1], c1, c2, 1.9]),color=(5/256,112/256,176/256),label="$u_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
# #ax_prions_2[0][1].plot(np.linspace(0,t[-1],100),np.array([(c1/c2)*(i+2)/(i+2) for i in range(100)]),"--",color=(0,0,0),label="$c_{1}/c_{2}$",linewidth=3)
# #ax_prions_2[1].plot(np.array([(1/(y0[0]*1.9-1))*np.log(1/vareps)*(i+2)/(i+2) for i in range(100)]),np.linspace(0,0.5,100),"--",color=(0,0,0),label="$\\tau=1/\\varepsilon$",linewidth=3)
# # Set a grid and define a legend
# ax_prions_2[1].grid()
# ax_prions_2[1].legend(loc='best',prop={"size":40})
# # Set the x-labels and y-labels
# ax_prions_2[1].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
# ax_prions_2[1].set_ylabel(ylabel="Prion conc.",fontsize=40)
# ax_prions_2[1].set_title(label="Healthy species",fontsize=50)    
# ax_prions_2[1].xaxis.set_tick_params(labelsize=35)
# ax_prions_2[1].yaxis.set_tick_params(labelsize=35)    
# #plt.show()
# f_prions_2.savefig('./Figures/HeMiTo_dynamics_2.png')


# # Non-linear conversion term
# plot_LaTeX_2D(t,u_2,"./Figures/prion_overall_dynamics/Input/healthy_eps_large.tex","color=blue_1,line width=1.0pt,","$u$")
# plot_LaTeX_2D(t,v_2,"./Figures/prion_overall_dynamics/Input/toxic_eps_large.tex","color=green_1,line width=1.0pt,","$v$")
# plot_LaTeX_2D(t2,u_He(t2, [y0[0], c1, c2])+vareps*u_Mi(t2, [y0[0], y0[1], c1, c2, f0]),"./Figures/prion_overall_dynamics/Input/healthy_eps_large.tex","color=blue_3,densely dashed,line width=1.0pt,","$u_{\\mathrm{He}}+\\varepsilon u_{\\mathrm{Mi}}$")
# plot_LaTeX_2D(t2,y0[1]*((t2+1)/(t2+1))+vareps*v_Mi(t2, [y0[0], y0[1], c1, c2, f0]),"./Figures/prion_overall_dynamics/Input/toxic_eps_large.tex","color=green_3,densely dashed,line width=1.0pt,","$v_{\\mathrm{He}}+\\varepsilon v_{\\mathrm{Mi}}$")


# # Linear conversion term
# plot_LaTeX_2D(t,u_1,"./Figures/prion_overall_dynamics_constant/Input/healthy_eps_small.tex","color=blue_1,line width=1.0pt,","$u$")
# plot_LaTeX_2D(t,v_1,"./Figures/prion_overall_dynamics_constant/Input/toxic_eps_small.tex","color=green_1,line width=1.0pt,","$v$")
# plot_LaTeX_2D(t2,u_He(t2, [y0[0], c1, c2])+vareps*u_Mi(t2, [y0[0], y0[1], c1, c2, 1]),"./Figures/prion_overall_dynamics_constant/Input/healthy_eps_small.tex","color=blue_3,densely dashed,line width=1.0pt,","$u_{\\mathrm{He}}+\\varepsilon u_{\\mathrm{Mi}}$")
# plot_LaTeX_2D(t2,y0[1]*((t2+1)/(t2+1))+vareps*v_Mi(t2, [y0[0], y0[1], c1, c2, 1]),"./Figures/prion_overall_dynamics_constant/Input/toxic_eps_small.tex","color=green_3,densely dashed,line width=1.0pt,","$v_{\\mathrm{He}}+\\varepsilon v_{\\mathrm{Mi}}$")





# # Parameters yeah
# p_high = [c1, 5.5*c2, vareps]
# print("High:\t\tc1/c2-u0\t=\t%0.3f"%(c1/(5.5*c2)-y0[0]))
# print("High c2:\t\tc2\t=\t%0.3f"%(5.5*c2))
# p_medium = [c1, 5.0*c2, vareps]
# print("Medium:\t\tc1/c2-u0\t=\t%0.3f"%(c1/(5.0*c2)-y0[0]))
# print("Medium c2:\t\tc2\t=\t%0.3f"%(5.0*c2))
# p_low = [c1, 4.5*c2, vareps]
# print("Low:\t\tc1/c2-u0\t=\t%0.3f"%(c1/(4.5*c2)-y0[0]))
# print("Low c2:\t\tc2\t=\t%0.3f"%(4.5*c2))
# # Define the age or time if you will
# t_exp = np.linspace(0,10,200)
# # Simulate patient with simple conversion term
# prions_exp_high = odeint(aB_Hill, y0, t_exp, args=(p_high,))
# prions_exp_medium = odeint(aB_Hill, y0, t_exp, args=(p_medium,))
# prions_exp_low = odeint(aB_Hill, y0, t_exp, args=(p_low,))
# # Extract each species for patient
# u_exp_high, v_exp_high = prions_exp_high.T
# u_exp_medium, v_exp_medium = prions_exp_medium.T
# u_exp_low, v_exp_low = prions_exp_low.T
# # Errors
# err_num_high = np.abs(v_exp_high-y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp))
# err_num_medium = np.abs(v_exp_medium-y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp))
# err_num_low = np.abs(v_exp_low-y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp))
# err_theo_high = vareps*y0[1]*(f0*(c1-5.5*c2*y0[0])-vareps*(y0[0]*f0-1)**2)*((t_exp**2)/2)
# err_theo_medium = vareps*y0[1]*(f0*(c1-5.0*c2*y0[0])-vareps*(y0[0]*f0-1)**2)*((t_exp**2)/2)
# err_theo_low = vareps*y0[1]*(f0*(c1-4.5*c2*y0[0])-vareps*(y0[0]*f0-1)**2)*((t_exp**2)/2)


# #Define the first figure
# f_prions_3, ax_prions_3 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# # Plot the simulated data and the underlying model
# ax_prions_3[0].plot(t_exp,v_exp_low,color=(35/256,139/256,69/256),label="$v(\\tau)$ low $c_{2}$",linewidth=3.0)
# ax_prions_3[0].plot(t_exp,v_exp_medium,color=(0/256,109/256,44/256),label="$v(\\tau)$ medium $c_{2}$",linewidth=3.0)
# ax_prions_3[0].plot(t_exp,v_exp_high,color=(0/256,68/256,27/256),label="$v(\\tau)$ high $c_{2}$",linewidth=3.0)
# # ax_prions_3[0].plot(t_exp,u_exp_high,color=(2/256,56/256,88/256),label="$u_{\\mathrm{high}}(\\tau)$",linewidth=3.0)
# # ax_prions_3[0].plot(t_exp,u_exp_medium,color=(4/256,90/256,141/256),label="$u_{\\mathrm{medium}}(\\tau)$",linewidth=3.0)
# # ax_prions_3[0].plot(t_exp,u_exp_low,color=(5/256,112/256,176/256),label="$u_{\\mathrm{low}}(\\tau)$",linewidth=3.0)
# ax_prions_3[0].plot(t_exp,y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp),color=(153/256,216/256,201/256),label="$v_{0}e^{\\lambda\\varepsilon\\tau}$",linewidth=3.0)
# # Set a grid and define a legend
# ax_prions_3[0].grid()
# ax_prions_3[0].legend(loc='best',prop={"size":40})
# # Set the x-labels and y-labels
# ax_prions_3[0].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
# ax_prions_3[0].set_ylabel(ylabel="Prion conc.",fontsize=40)
# ax_prions_3[0].set_title(label="Toxic species",fontsize=50)    
# ax_prions_3[0].xaxis.set_tick_params(labelsize=35)
# ax_prions_3[0].yaxis.set_tick_params(labelsize=35)
# # Plot the simulated data and the underlying model
# ax_prions_3[0].plot(t_exp,v_exp_low,color=(35/256,139/256,69/256),label="$v_{\\mathrm{low}}(\\tau)$",linewidth=3.0)
# #========================================================================
# # Second plot
# #========================================================================
# ax_prions_3[1].plot(t_exp,100*np.abs(err_theo_high)/(np.abs(err_theo_high)+np.abs(err_num_high-err_theo_high)),color=(35/256,139/256,69/256),label="$R_{\\mathrm{HSS}}(\\tau)$ high $c_{2}$",linewidth=3.0)
# ax_prions_3[1].plot(t_exp,100*np.abs(err_theo_medium)/(np.abs(err_theo_medium)+np.abs(err_num_medium-err_theo_medium)),color=(0/256,109/256,44/256),label="$R_{\\mathrm{HSS}}(\\tau)$ medium $c_{2}$",linewidth=3.0)
# ax_prions_3[1].plot(t_exp,100*np.abs(err_theo_low)/(np.abs(err_theo_low)+np.abs(err_num_low-err_theo_low)),color=(0/256,68/256,27/256),label="$R_{\\mathrm{HSS}}(\\tau)$ low $c_{2}$",linewidth=3.0)
# #ax_prions_3[1].plot(t_exp,err_theo_high,color=(77/256,0/256,75/256),label="$E_{\\mathrm{HSS}}(\\tau)$ high $c_{2}$",linewidth=3.0)
# #ax_prions_3[1].plot(t_exp,err_theo_medium,color=(129/256,15/256,124/256),label="$E_{\\mathrm{HSS}}(\\tau)$ medium $c_{2}$",linewidth=3.0)
# #ax_prions_3[1].plot(t_exp,err_theo_low,color=(136/256,65/256,157/256),label="$E_{\\mathrm{HSS}}(\\tau)$ low $c_{2}$",linewidth=3.0)
# ax_prions_3[1].grid()
# ax_prions_3[1].legend(loc='best',prop={"size":40})
# # Set the x-labels and y-labels
# ax_prions_3[1].set_xlabel(xlabel="$\\tau$",fontsize=40)
# ax_prions_3[1].set_ylabel(ylabel="$E(\\tau)$",fontsize=40)
# ax_prions_3[1].set_title(label="Error in terms of degradation",fontsize=50)    
# ax_prions_3[1].xaxis.set_tick_params(labelsize=35)
# ax_prions_3[1].yaxis.set_tick_params(labelsize=35)
# plt.show()
# f_prions_3.savefig('./Figures/HeMiTo_dynamics_3.png')


# # Exponential growth
# plot_LaTeX_2D(t_exp,v_exp_low,"./Figures/prion_exponential/Input/v_exp.tex","color=green_3,line width=1.0pt,loosely dotted,mark=triangle*,mark repeat=20","$\left.v(\\tau)\\right|_{c_{2}=1.35}$")
# plot_LaTeX_2D(t_exp,v_exp_medium,"./Figures/prion_exponential/Input/v_exp.tex","color=green_2,line width=1.0pt,densely dashdotted,mark=square*,mark repeat=20","$\left.v(\\tau)\\right|_{c_{2}=1.50}$")
# plot_LaTeX_2D(t_exp,v_exp_high,"./Figures/prion_exponential/Input/v_exp.tex","color=green_1,line width=1.0pt, densely dotted, mark=square*,mark repeat=20","$\left.v(\\tau)\\right|_{c_{2}=1.65}$")
# plot_LaTeX_2D(t_exp,y0[1]*np.exp((y0[0]*f0-1)*vareps*t_exp),"./Figures/prion_exponential/Input/v_exp.tex","color=green_4,line width=1.0pt,densely dashed, mark=otimes*,mark repeat=20","$v_{0}e^{\\lambda\\varepsilon\\tau}$")


# #  Errors, numerical and theoretical
# # Low c2
# #plot_LaTeX_2D(t_exp,err_num_low,"./Figures/prion_exponential/Input/err_exp.tex","color=green_3,line width=1.0pt,loosely dotted,mark=triangle*,mark repeat=20","$\left.|v(\\tau)-v_{0}e^{\\lambda\\varepsilon\\tau}|\\right|_{c_{2}=1.35}$")
# plot_LaTeX_2D(t_exp,100*np.abs(err_theo_low)/(np.abs(err_theo_low)+np.abs(err_num_low-err_theo_low)),"./Figures/prion_exponential/Input/err_exp.tex","color=magenta_3,line width=1.0pt,loosely dotted,mark=square*,mark repeat=20","$\left.R_{\\mathrm{HSS}}(\\tau)\\right|_{c_{2}=1.35}$")
# # Medium c2
# #plot_LaTeX_2D(t_exp,err_num_medium,"./Figures/prion_exponential/Input/err_exp.tex","color=green_2,line width=1.0pt,densely dashdotted,mark=square*,mark repeat=20","$\left.|v(\\tau)-v_{0}e^{\\lambda\\varepsilon\\tau}|\\right|_{c_{2}=1.50}$")
# plot_LaTeX_2D(t_exp,100*np.abs(err_theo_medium)/(np.abs(err_theo_medium)+np.abs(err_num_medium-err_theo_medium)),"./Figures/prion_exponential/Input/err_exp.tex","color=magenta_2,line width=1.0pt,densely dashdotted,mark=halfcircle*,mark repeat=20","$\left.R_{\\mathrm{HSS}}(\\tau)\\right|_{c_{2}=1.50}$")
# # High c2
# #plot_LaTeX_2D(t_exp,err_num_high,"./Figures/prion_exponential/Input/err_exp.tex","color=green_1,line width=1.0pt, densely dotted, mark=halfsquare right*,mark repeat=20","$\left.|v(\\tau)-v_{0}e^{\\lambda\\varepsilon\\tau}|\\right|_{c_{2}=1.65}$")
# plot_LaTeX_2D(t_exp,100*np.abs(err_theo_high)/(np.abs(err_theo_high)+np.abs(err_num_high-err_theo_high)),"./Figures/prion_exponential/Input/err_exp.tex","color=magenta_1,line width=1.0pt, densely dotted, mark=halfdiamond*,mark repeat=20","$\left.R_{\\mathrm{HSS}}(\\tau)\\right|_{c_{2}=1.65}$")



