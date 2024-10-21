import numpy as np
import matplotlib.pyplot as plt
import nuclei_func as nf
from scipy.interpolate import make_interp_spline

h_from_c = [42.5, 127.5, 212.5, 297.5, 382.5, 467.5, 552.5, 637.5, 722.5, 807.5]
peak_factors = [0.572, 0.737, 0.868, 0.958, 1, 0.983, 0.912, 0.802, 0.658, 0.498]
A = 238
sigma_scattering = nf.macro(2.5, 10.5, 238)
sigma_trans = sigma_scattering * (1 + 2 / (3 * A))
extrapolated_length = 0.714 / sigma_trans
print(extrapolated_length)

def power_profile(x):
    length_h_from_c = len(h_from_c)

    # Define a vector that returns the start and end points of each interval
    length_interval = 1 + length_h_from_c  # One extra for the initial 0, and one for the final 850
    interval = np.zeros(length_interval)

    # Set the first and last values explicitly
    interval[0] = 0
    interval[-1] = 850

    # Calculate midpoints for each interval
    for i in range(1, length_h_from_c):
        interval[i] = (h_from_c[i - 1] + h_from_c[i]) / 2

    # Define a function that has a cosine in each interval with the peak at the midpoint and varying amplitude
    value = 0
    for i in range(1, len(interval)):
        start = interval[i - 1]
        end = interval[i]
        midpoint = (start + end) / 2
        amplitude = peak_factors[i - 1]
        adjusted_start = start - extrapolated_length
        adjusted_end = end + extrapolated_length
        frequency = np.pi / (adjusted_end - adjusted_start)
        if adjusted_start <= x < adjusted_end:
            value += amplitude * np.cos(frequency * (x - midpoint))
    return [value, interval]

# Generate x values and corresponding y values for the plot
x_values = np.linspace(0, 850, 1000)
y_values = [power_profile(x)[0] for x in x_values]

# Use spline interpolation to smooth the curve without altering the peaks
interval = power_profile(x_values[0])[1]
x_peaks = [(interval[i - 1] + interval[i]) / 2 for i in range(1, len(interval))]
y_peaks = peak_factors
spl = make_interp_spline(x_values, y_values, k=3)
y_values_smooth = spl(x_values)

# Calculate maximum and minimum points
max_points_x = x_peaks
max_points_y = peak_factors
min_points_x = interval
min_points_y = [0] * len(interval)

# Plotting the graph
plt.figure(figsize=(10, 6))

# Plot dotted lines at interval start and end
for point in interval:
    plt.axvline(x=point, color='r', linestyle='--', linewidth=1)

# Plot dotted lines at midpoints
for point in h_from_c:
    plt.axvline(x=point, color='b', linestyle='-.', linewidth=1)

# Plot the smoothed cosine function
plt.plot(x_values, y_values_smooth, color='g', label='Smoothed Cosine Function in Each Interval with Peak Factors')

# Plot black dashed line connecting the maximum points
plt.plot(max_points_x, max_points_y, color='k', linestyle='--', linewidth=1.5, label='Max Points Connection')

# Plot black dashed line connecting the minimum points
plt.plot(min_points_x, min_points_y, color='k', linestyle='--', linewidth=1.5, label='Min Points Connection')

plt.xlabel('Position')
plt.ylabel('Value')
plt.title('Interval Start, End, Midpoints, and Smoothed Cosine Function with Peak Factors')
plt.grid(True)
plt.legend()
plt.show()