# 2D plot of the temperature profile

R = fuel_outer_diameter/2
R_void = 0.511873*10**-3
R_new, R_start, T_hot = cold_to_hot_fuel(Fuel_Proprieties,Geometrical_Data,vars,h_vals)
R_new_value = R_new[0]
print(R_new_value)
r_fuel_vector = np.linspace(R_void,R_new,1000)

r_gap_fuel = R_new_value
r_end = R_void

r_fuel = np.linspace(r_gap_fuel, r_end, 25)



temp_plot_bottom = [f.get_temperature_at_point(0, r, vars.T_map) for r in r_plot]
temp_plot_center = [f.get_temperature_at_point(0.425, r, vars.T_map) for r in r_plot]
temp_plot_top = [f.get_temperature_at_point(0.850, r, vars.T_map) for r in r_plot]

# Create plot
plt.plot(r_plot * 1e3, temp_plot_bottom, label='Bottom', marker='o')
plt.plot(r_plot * 1e3, temp_plot_center, label='Center', marker='o')
plt.plot(r_plot * 1e3, temp_plot_top, label='Top', marker='o')

r_0 = Geometrical_Data.fuel_inner_diameter/2 * 1e3
r_1 = Geometrical_Data.fuel_outer_diameter/2 * 1e3
r_2 = Geometrical_Data.cladding_outer_diameter/2 * 1e3 - Geometrical_Data.thickness_cladding * 1e3
r_3 = Geometrical_Data.cladding_outer_diameter/2 * 1e3
r_4 = r_plot[0]*1e3

# Add shading to different regions
colors = ['#00008B', '#0000CD', '#4169E1', '#6495ED', '#87CEEB']
plt.axvspan(r_3, r_4, color=colors[1], alpha=0.3, label='Coolant Region')
plt.axvspan(r_2, r_3, color=colors[2], alpha=0.3, label='Cladding Region')
plt.axvspan(r_1, r_2, color=colors[3], alpha=0.3, label='Gap Region')
plt.axvspan(r_0, r_1, color=colors[4], alpha=0.3, label='Fuel Region')

# Set title and axis labels
plt.title('Temperature Profile')
plt.xlabel('Radius [mm]')
plt.xlim(r_0, r_4)
plt.ylabel('Temperature [K]')

# Add legend to the plot
plt.legend()

# Put the legend out of the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show the figure
plt.show()