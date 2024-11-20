import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

def calcification_model(external_calcium_amount,bound_calmodulin_amount):
    # Calcium constants 
    external_cal_conc = external_calcium_amount # in mM (grown in marine water with 10mM Ca conc)
    cellular_cal_conc = 1e-4 # in mM (has 100 nM - maintained in this range as higher value is cytotoxic)
    weight_of_calcite = 100 #in gram/mol
    cocco_weight= 2 #in picogram per coccolith (average weight)
    calcite_percent= 0.86 #percentage of calcite in coccolith
    calmodulin_equilibrium_constant = 1e9
    pH_of_CV =  (math.log10(external_cal_conc/cellular_cal_conc)+11) / 2 #calcification modulates inhibition by changing pH of CV (range between 7 and 8)
    if pH_of_CV < 7:
        pH_of_CV = 7
    if pH_of_CV > 8:
        pH_of_CV = 8
    total_calcification =0 #in pmol per cell
    total_C_fixed =0 #in pg per cell
    total_coccolith_produced =0 #per cell over its lifespan of 72 hours 
    storage_amount =0 # storage gets filled if uptake greater than deposition capacity
    # (Here storage acts as just a placeholder, such that when uptake exceeds deposition capacity feedback takes place to adjust
    # uptake, though in this model instead of dynamically changing uptake rate depending on this internal feedback, it is taken 
    # based on external, internal calcium conc. and free calmodulin amount only. But the feedback gets adjusted since that amount 
    # doesnt contribute to final coccolith formation so instead of next time cycle getting affected effect nullifies in current 
    # time period)

    bound_calmodulin = bound_calmodulin_amount # in M concentration (see notes for explanation)
    free_calmodulin = bound_calmodulin / (calmodulin_equilibrium_constant*cellular_cal_conc*1e-3) # factor is multiplied for unit matching 
    uptake_rate = (2 * (math.log10(external_cal_conc/cellular_cal_conc)-3) * math.exp(1e7 * free_calmodulin)) / 24# in pg per cell per hour
    cal_cellular_need = 0.5 / 24 # in pg per cell per hour
    amount_of_cal_calcification = uptake_rate-cal_cellular_need # in pg per cell per hour
    calcification_uptake_hourly = amount_of_cal_calcification/40 # in pmol per cell per hour

    # Constants given
    change_gpa_assembly_constant = 0.1
    nucleators_per_GPA = 1
    calcite_deposition_constant = 1
    degradation_inhibitor_constant = 0.015
    nucleator_reduction_constant = 1e-12
    a = 10
    b = 0.05
    gpa_assembly_rate_0 = 1

    # Initial conditions
    initial_gpa = 0
    initial_inhibitor = 0.5
    intial_nucleator = 0
    initial_deposition_capacity = 0

    # Function that returns changes in parameters - differential equation 
    def model(z, t):
        gpa,inhibitor, nucleator, deposition_capacity = z
        ddeposition_capacity_dt = calcite_deposition_constant * (b / (b + inhibitor)**a) * nucleator
        dinhibitor_dt = - degradation_inhibitor_constant * gpa* inhibitor * (pH_of_CV-6.9999)
        dgpa_dt = gpa_assembly_rate_0 - change_gpa_assembly_constant * gpa
        dnucleator_dt = nucleators_per_GPA * dgpa_dt - nucleator_reduction_constant * ddeposition_capacity_dt * nucleator
        return [dgpa_dt, dinhibitor_dt, dnucleator_dt, ddeposition_capacity_dt]

    # Time points
    t = np.linspace(0, 72, 72*10)  # from 0 to 72 hours, 720 points

    # Initial conditions vector
    initial_conditions = [initial_gpa, initial_inhibitor, intial_nucleator, initial_deposition_capacity]

    # Solve ODE
    z = odeint(model, initial_conditions, t)

    # Extracting results
    gpa = z[:, 0]
    inhibitor = z[:, 1]
    nucleator = z[:, 2]
    deposition_capacity = z[:, 3]

    selected_indices = [i for i in range(0, 720, 10)]
    gpa_hourly = [gpa[i] for i in selected_indices]
    inhibitor_hourly = [inhibitor[i] for i in selected_indices]
    nucleator_hourly = [nucleator[i] for i in selected_indices]
    deposition_capacity_hourly = [deposition_capacity[i] for i in selected_indices] #in molecules per micrometer^3

    volume_of_deposition = 4 * np.pi * (pow(3.5,3)-pow(3,3)) / 3 # cell radius 3 micrometer and consider 5 layers of 0.1 micrometer thickness each layer
    deposition_capacity_final = np.array(deposition_capacity_hourly) * (volume_of_deposition * 1e12 / (6.022*1e23)) # in pmol per hour per cell

    #Replicating the process - uptake indicates availability, CV pH dictates inhibition to deposition process which is in turn dependant 
    # concentration gradient and deposition capability takes into account final deposition with storage vesicle taking care when uptake
    # is higher than deposition capability
    for i in range(72): #over 72 hours lifespan
        if deposition_capacity_final[i] < calcification_uptake_hourly :
            total_calcification+= deposition_capacity_final[i]
            residue = calcification_uptake_hourly-deposition_capacity_final[i]
            storage_amount = storage_amount+residue 
        else:
            total_calcification+= calcification_uptake_hourly

    total_C_fixed = total_calcification*12 # in pg per cell per hour
    total_coccolith_produced = math.floor(total_calcification*weight_of_calcite/(cocco_weight*calcite_percent))
    if total_calcification<=0 :
        total_calcification=0
        total_C_fixed=0
        total_coccolith_produced=0

    return [total_calcification,total_C_fixed,total_coccolith_produced]

    '''
    # for understanding that uptake is limiting than deposition
    calci_hourly = np.full(72, calcification_uptake_hourly)
    production_overall = np.minimum(calci_hourly, deposition_capacity_final)

    cumulative_calci_hourly = np.array([])
    cumulative_sum = 0
    for value in calci_hourly:
        cumulative_sum += value
        cumulative_calci_hourly = np.append(cumulative_calci_hourly, cumulative_sum)

    cumulative_deposition = np.array([])
    cumulative_sum = 0
    for value in deposition_capacity_final:
        cumulative_sum += value
        cumulative_deposition = np.append(cumulative_deposition, cumulative_sum)

    cumulative_production = np.array([])
    cumulative_sum = 0
    for value in production_overall:
        cumulative_sum += value
        cumulative_production = np.append(cumulative_production, cumulative_sum)
    lifespan = np.linspace(0, 72, 72) 
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(lifespan[:25], cumulative_deposition[:25], label='Cumulative Deposition capacity', color='blue')         
    plt.plot(lifespan[:25], cumulative_calci_hourly[:25], label='Cumulative Calcification based on uptake', color='grey')   
    plt.plot(lifespan[:25], cumulative_production[:25], label='Cumulative Final amount deposited', color='red')  

    plt.xlabel('Time (in hrs)')
    plt.ylabel('Amount of Calcium Carbonate (in pmol per hour per cell)')
    plt.title('To show that uptake is limiting over deposition using cumulative plot(zoomed out)')
    plt.legend()
    plt.grid()
    plt.savefig('cumulative_compare_hourly_zoomed.png')  
    '''

if __name__ == "__main__":
    external_calci = []
    external_cocco = []
    external_C_fixed = []
    external_val = []

    bound_calci = []
    bound_cocco = []
    bound_C_fixed = []
    bound_val = []
    
    for i in range(1,101):
        external_val.append(i*0.1) #external calcium conc. is varied from 0.1 mM to 10 mM (in sea water it is 10 mM)
        returned_arr = calcification_model(i*0.1,1e-5) #calmodulin amount is as expressed in normal cell 
        external_calci.append(returned_arr[0])
        external_C_fixed.append(returned_arr[1])
        external_cocco.append(returned_arr[2])
    
    for i in range(1,101):
        bound_val.append((1e-5+i*1e-7)*1e3) #as bound val is increased (actually expression is increasing by same proportion)
        # The free CAM also increasing by same proportion for maintining equilibrium value thus increasing activation
        returned_arr = calcification_model(10,1e-5+i*1e-7) #calmodulin amount is as expressed in normal cell 
        bound_calci.append(returned_arr[0])
        bound_C_fixed.append(returned_arr[1])
        bound_cocco.append(returned_arr[2])
    
    temp_arr = calcification_model(10,1e-5)
    initial_val =temp_arr[0]

    temp_arr = calcification_model(10*1.1,1e-5)
    val = temp_arr[0]
    percent_increase = (val-initial_val)*100/initial_val
    print("Increase in 10 percent of external calcium above normal sea level leads to {:.2f} percent increase in calcification".format(percent_increase))

    temp_arr = calcification_model(10,(1e-5) *1.1)
    val = temp_arr[0]
    percent_increase = (val-initial_val)*100/initial_val
    print("Increase in 10 percent of calmodulin expression above normal levels leads to {:.2f} percent increase in calcification".format(percent_increase))

    fig, axs = plt.subplots(1, 3, figsize=(15, 10))  

    axs[0].plot(external_val, external_calci, label='Total Calcification')
    axs[0].set_title('Total Calcification wrt changing external calcium conc.')
    axs[0].set_xlabel('External Calcium Conc. (in mM)')
    axs[0].set_ylabel('Calcification (in pmol per cell)')
    axs[0].legend()

    axs[1].plot(external_val, external_C_fixed, label='Carbon Fixed')
    axs[1].set_title('Carbon fixed wrt changing external calcium conc.')
    axs[1].set_xlabel('External Calcium Conc. (in mM)')
    axs[1].set_ylabel('Amount of Carbon fixed (in pg per cell)')
    axs[1].legend()

    axs[2].plot(external_val, external_cocco, label='Coccoliths produced')
    axs[2].set_title('Number of coccoliths produced wrt changing external calcium conc.')
    axs[2].set_xlabel('External Calcium Conc. (in mM)')
    axs[2].set_ylabel('Coccoliths produced')
    axs[2].legend()

    plt.tight_layout()
    plt.savefig('calci_external_graph.png')

    fig, axs = plt.subplots(1, 3, figsize=(15, 10))  

    axs[0].plot(bound_val, bound_calci, label='Total Calcification')
    axs[0].set_title('Total Calcification wrt changing calmodulin expression')
    axs[0].set_xlabel('Bound Calmodulin (in mM) ')
    axs[0].set_ylabel('Calcification (in pmol per cell)')
    axs[0].legend()

    axs[1].plot(bound_val, bound_C_fixed, label='Carbon Fixed')
    axs[1].set_title('Carbon fixed wrt changing calmodulin expression')
    axs[1].set_xlabel('Bound Calmodulin (in mM) ')
    axs[1].set_ylabel('Amount of Carbon fixed (in pg per cell)')
    axs[1].legend()

    axs[2].plot(bound_val, bound_cocco, label='Coccolith produced')
    axs[2].set_title('Number of coccoliths produced wrt changing calmodulin expression')
    axs[2].set_xlabel('Bound Calmodulin (in mM) ')
    axs[2].set_ylabel('Coccoliths produced')
    axs[2].legend()

    plt.tight_layout()
    plt.savefig('calci_calmodulin_graph.png')

    #Normalization to remove amount factor and see percentwise increase in calcification (and thus carbon fixed and coccolith produced)
    # because of changing external calcium concentration and calmodulin expressed above normal levels of 10mM and 1e-5 M 
    
    percent_change = []
    percent_change_external = []
    percent_change_bound = []

    for i in range(0,101,5):
        percent_change.append(i)
        temp_arr = calcification_model(10*(1+0.01*i),1e-5)
        percent_change_external.append((temp_arr[0]-initial_val)* 100 / initial_val)
        temp_arr = calcification_model(10,1e-5*(1+0.01*i))
        percent_change_bound.append((temp_arr[0]-initial_val)* 100 /initial_val)

    plt.figure(figsize=(10, 6))
    plt.plot(percent_change, percent_change_external, label='Percent change in calcification wrt to percent change in external calcium conc.', color='blue')         
    plt.plot(percent_change, percent_change_bound, label='Percent change in calcification wrt to percent change in calmodulin', color='red')  

    plt.xlabel('Percent Change in variable')
    plt.ylabel('Percent Change in calcification')
    plt.title('To show modulating calmodulin is effective')
    plt.legend()
    plt.grid()
    plt.savefig('external_vs_CAM.png')  



    

    




