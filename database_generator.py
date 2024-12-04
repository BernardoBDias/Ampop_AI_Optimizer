import os
from PyLTSpice import SimCommander
import ltspice
import numpy as np
import matplotlib.pyplot as plt
import csv

n_amostras = 100

# Criação e ediação da netlist:
my_path = os.path.dirname(os.path.realpath(__file__))
# print('rel path: ', my_path)
LTcircuit_OL = SimCommander(my_path+"\\Ampop_sim_malha_aberta.asc")
LTcircuit_SR = SimCommander(my_path+"\\Ampop_sim_SR.asc")

filename = "database_"+str(n_amostras)+"_samples.csv"
headers = ["Wd", "Wn", "Wout", "N", "DC Gain", "Freq Corte", "SR"]

with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(headers)

#Parâmetros a serem simulados:
W_d = np.random.randint(low=1,high=50,size=n_amostras)*0.2e-6 
W_n = np.random.randint(low=1,high=50,size=n_amostras)*0.2e-6 
W_out = np.random.randint(low=1,high=50,size=n_amostras)*0.2e-6 
M_n = np.random.randint(low=1,high=20,size=n_amostras)

for trial in range(n_amostras):
    LTcircuit_OL.set_parameters(Wd = W_d[trial])
    LTcircuit_OL.set_parameters(Wn = W_n[trial])
    LTcircuit_OL.set_parameters(Wout = W_out[trial])
    LTcircuit_OL.set_parameters(N = M_n[trial])
    LTcircuit_SR.set_parameters(Wd = W_d[trial])
    LTcircuit_SR.set_parameters(Wn = W_n[trial])
    LTcircuit_SR.set_parameters(Wout = W_out[trial])
    LTcircuit_SR.set_parameters(N = M_n[trial])
    LTcircuit_OL.run()
    LTcircuit_OL.wait_completion()
    LTcircuit_SR.run()
    LTcircuit_SR.wait_completion()

    # Análise de dados:
    # circuit_file = "LP_test_1.raw"
    # absolute_path = my_path+'\\'+circuit_file
    # print('abs path: ', absolute_path)
    # circuit_data = ltspice.Ltspice(os.path.dirname(__file__)+'\\Ampop_sim_malha_aberta_1.raw')
    circuit_data = ltspice.Ltspice(os.path.dirname(__file__)+'\\Ampop_sim_malha_aberta_'+str(trial+1)+'.raw')
    circuit_data.parse()
    # print('Parse and raw done!')

    freq = circuit_data.get_frequency()
    V_in = circuit_data.get_data('V(V2)')
    V_out = circuit_data.get_data('V(Vo)')

    # print('Vout len: ', len(np.real(V_out)))
    # print('Freq len: ', len(freq))

    n_points = len(np.real(V_out))
    dc_value = abs(V_out[0])/abs(V_in[0])
    for i in range(n_points):
        if abs(V_out[i])/abs(V_in[i]) <= dc_value/np.sqrt(2):
            print('-3dB = ', abs(V_out[i])/abs(V_in[i]))
            band_pass = freq[i]
            break
    
    circuit_data = ltspice.Ltspice(os.path.dirname(__file__)+'\\Ampop_sim_SR_'+str(trial+1)+'.raw')
    # circuit_data = ltspice.Ltspice(os.path.dirname(__file__)+'\\Ampop_sim_SR_1.raw')
    circuit_data.parse()

    V_out = circuit_data.get_data('V(Vo)')
    sim_time_data = circuit_data.get_time()
    deriv = np.gradient(V_out,sim_time_data)
    slew_rate = max(abs(deriv))

    trial_data = [W_d[trial],W_n[trial],W_out[trial],M_n[trial],dc_value,band_pass, slew_rate]
    with open(filename, mode='a', newline='') as file:  # Abre o arquivo em modo de adição
        writer = csv.writer(file)
        writer.writerow(trial_data)

    # print('Frequencia de corte: ', band_pass)
    # print('Ganho DC: ', dc_value)
    # print('SR: ', slew_rate)
    if n_amostras/trial == 10:
        print('10%')
    if n_amostras/trial == 5:
        print('20%')
    if n_amostras/trial == 2:
        print('50%')
    if n_amostras/trial == 1.25:
        print('80%')

    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_malha_aberta_'+str(trial+1)+'.raw')
    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_malha_aberta_'+str(trial+1)+'.op.raw')
    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_malha_aberta_'+str(trial+1)+'.net')
    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_malha_aberta_'+str(trial+1)+'.log')
    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_SR_'+str(trial+1)+'.raw')
    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_SR_'+str(trial+1)+'.op.raw')
    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_SR_'+str(trial+1)+'.net')
    os.remove(os.path.dirname(__file__)+'\\Ampop_sim_SR_'+str(trial+1)+'.log')

print('Done!')
# resistencias = [500, 1e3, 2e3]

# for r in resistencias:
#     LTcircuit.set_parameters(R_var = r)
#     print('R = ', r)
#     LTcircuit.run()
#     LTcircuit.wait_completion()

#     # Análise de dados:
#     # circuit_file = "LP_test_1.raw"
#     # absolute_path = my_path+'\\'+circuit_file
#     # print('abs path: ', absolute_path)
#     circuit_data = ltspice.Ltspice(os.path.dirname(__file__)+'\\LP_test_'+str(resistencias.index(r)+1)+'.raw')
#     circuit_data.parse()
#     print('Parse and raw done!')

#     freq = circuit_data.get_frequency()
#     V_in = circuit_data.get_data('V(Vin)')
#     V_out = circuit_data.get_data('V(Vout)')

#     print('Vout len: ', len(np.real(V_out)))
#     print('Freq len: ', len(freq))

#     n_points = len(np.real(V_out))
#     dc_value = abs(V_out[0])/abs(V_in[0])
#     for i in range(n_points):
#         if abs(V_out[i])/abs(V_in[i]) <= dc_value/np.sqrt(2):
#             print('-3dB = ', abs(V_out[i])/abs(V_in[i]))
#             band_pass = freq[i]
#             break

#     print('Frequencia de corte: ', band_pass)

    # plt.figure('AC')
    # plt.subplot(1,2,1)
    # plt.semilogx(freq, np.real(20*np.log10((V_out/V_in))))
    # plt.subplot(1,2,2)
    # plt.semilogx(freq, np.rad2deg(20*np.angle(V_out/V_in)))
    # plt.show()


