import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker



#Initialise the size of the box
size = np.array([3, 11, 3])#[x cm, y cm, z cm]
step = 0.01#cm
grid_size = np.ceil(size/step).astype('int')

#Calculate the centre of the box
origin = [grid_size[0]/2-0.5, grid_size[1]/2-0.5, grid_size[2]/2-0.5]#[x, y, z]

#sdfsdfsdfsdfsdfsdfsdf
x_distance = np.arange(grid_size[0])*step
y_distance = np.arange(grid_size[1])*step
z_distance = np.arange(grid_size[2])*step

print(grid_size)
print(x_distance)
print(x_distance-0.5-grid_size[0]/2)
print(y_distance-y_distance[-1]/2)
print(z_distance-z_distance[-1]/2)







M = np.array([[0], [0.000000015], [0]])




def calculate_angles(xoffset, yoffset, mode):

	if(mode == 'degrees'):
		#Convert degrees to radians
		xoffset = np.radians(xoffset)
		yoffset = np.radians(yoffset)

	#Create rotation matrices
	x_rotation_matrix = np.array([[1, 0, 0], [0, np.cos(xoffset), -np.sin(xoffset)], [0, np.sin(xoffset), np.cos(xoffset)]])
	y_rotation_matrix = np.array([[np.cos(yoffset), 0, np.sin(yoffset)], [0, 1, 0], [-np.sin(yoffset), 0, np.cos(yoffset)]])
	
    #Calculate the dipole vector
	dipole_vertical_vector = np.dot(y_rotation_matrix, np.dot(x_rotation_matrix, np.array([[0], [size[1]], [0]])))
	dipole_horizontal_vector = np.dot(np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]]), dipole_vertical_vector)

	coordinate_grid = np.zeros((grid_size[0], grid_size[1], grid_size[2], 3))
	cos_theta = np.zeros((grid_size[0], grid_size[1], grid_size[2]))
	cos_phi = np.zeros((grid_size[0], grid_size[1], grid_size[2]))

	for i in range(0,grid_size[0]):
		for j in range(0,grid_size[1]):
			#Create a 3D grid with the coordinates of each point stored in a 1D array at each point
			coordinate_grid[i][j] = np.array((grid_size[0]/2-0.5-i, grid_size[1]/2-0.5-j, 0))*step
			coordinate_grid[i][j][:,2] = (grid_size[2]/2-0.5-np.arange(grid_size[2]))*step
			for k in range(0,grid_size[2]):
				#Calculate the cosine of the angle between the dipole vector and each point
				cos_theta[i][j][k] = np.dot(dipole_vertical_vector.T, coordinate_grid[i][j][k])/(np.linalg.norm(dipole_vertical_vector)*np.linalg.norm(coordinate_grid[i][j][k]))
				cos_phi[i][j][k] = np.dot(dipole_horizontal_vector.T, coordinate_grid[i][j][k])/(np.linalg.norm(dipole_horizontal_vector)*np.linalg.norm(coordinate_grid[i][j][k]))
	#g = np.full((grid_size[0], grid_size[1], grid_size[2], 3), dipole_vertical_vector)
	return np.arccos(cos_theta), np.arccos(cos_phi)











#Create a 3D array with the euclidean distance from the origin stored at each point in metres
x, y, z = np.ogrid[0:grid_size[0], 0:grid_size[1],0: grid_size[2]]
r = (((x-origin[0])**2 + (y-origin[1])**2 + (z-origin[2])**2)**0.5)*step
r_cubed = (((x-origin[0])**2 + (y-origin[1])**2 + (z-origin[2])**2)**1.5)*step/1000000
print(r)




############################################################################################################################################################












#Create a 3D arrays of the angle between the vertical and horizontal vectors and each point stored at each point
theta, phi = calculate_angles(90, 0, 'degrees')


#Calculate the magnetic field strength at each point in the grid
Bx = 3e-7*np.linalg.norm(M)*np.sin(theta)*np.cos(theta)*np.cos(phi)/r_cubed
By = 3e-7*np.linalg.norm(M)*np.sin(theta)*np.cos(theta)*np.sin(phi)/r_cubed
Bz = 1e-7*np.linalg.norm(M)*(3*(np.cos(theta))**2-1)/r_cubed

#np.set_printoptions(threshold=np.inf)
#print(np.round(r_cubed**(1/3),5))	points+=1
#print(np.round(theta*180/np.pi,5))
#print(np.round(phi*180/np.pi,5))

ring_radius = 5
ring_seperation = 2

circle_mask = np.full((grid_size[0], grid_size[2]), 1)
area_element = grid_size[0] * grid_size[2]
for i in range(grid_size[0]):
	for j in range(grid_size[2]):
		if((i-origin[0])**2 + (j-origin[2])**2 > ring_radius**2):
			circle_mask[i][j] = 0
			area_element -= 1

flux = np.zeros(grid_size[1] + (2*int(ring_seperation/step))+4)
z_position = (np.arange(grid_size[1] + (2*int(ring_seperation/step))+4) - grid_size[1]/2 - (int(ring_seperation/step))+2)/10

for i in range(grid_size[1] + (2*int(ring_seperation/step))+4):
	if(i >= 0 and i < grid_size[1]):
		Bx[:, i, :] *= circle_mask
		cos_theta = np.zeros((grid_size[0], grid_size[2]))

		for j in range(grid_size[0]):
			for k in range(grid_size[2]):
				B_vector = np.array((Bx[j][i][k], By[j][i][k], Bz[j][i][k]))
				cos_theta[j][k] = np.dot(B_vector, np.array([[0], [size[1]], [0]]))/(np.linalg.norm(B_vector)*np.linalg.norm(np.array([[0], [size[1]], [0]])))
				if(np.isnan(cos_theta[j][k])):
					cos_theta[j][k] = 0

		flux[i] += ((np.sum(Bx[:, i, :])**2 + np.sum(By[:, i, :])**2 + np.sum(Bz[:, i, :])**2)**0.5) * (np.pi/area_element)*ring_radius**2 * abs(np.sum(cos_theta))

	if(int(i - (ring_seperation/step)+2) >= 0 and int(i - (ring_seperation/step)+2) < grid_size[1]):
		cos_theta = np.zeros((grid_size[0], grid_size[2]))

		for j in range(grid_size[0]):
			for k in range(grid_size[2]):
				B_vector = np.array((Bx[j][int(i - (ring_seperation/step)+2)][k], By[j][int(i - (ring_seperation/step)+2)][k], Bz[j][int(i - (ring_seperation/step)+2)][k]))
				cos_theta[j][k] = np.dot(B_vector, np.array([[0], [size[1]], [0]]))/(np.linalg.norm(B_vector)*np.linalg.norm(np.array([[0], [size[1]], [0]])))
				if(np.isnan(cos_theta[j][k])):
					cos_theta[j][k] = 0

		flux[i] -= 2*((np.sum(Bx[:, int(i - (ring_seperation/step)+2), :])**2 + np.sum(By[:, int(i - (ring_seperation/step)+2), :])**2 + np.sum(Bz[:, int(i - (ring_seperation/step)+2), :])**2)**0.5) * (np.pi/area_element)*ring_radius**2 * abs(np.sum(cos_theta))

	if(i - (2*(ring_seperation/step))+4 >= 0 and i - (2*(ring_seperation/step))+4 < grid_size[1]):
		cos_theta = np.zeros((grid_size[0], grid_size[2]))

		for j in range(grid_size[0]):
			for k in range(grid_size[2]):
				B_vector = np.array((Bx[j][i - (2*int(ring_seperation/step))+4][k], By[j][i - (2*int(ring_seperation/step))+4][k], Bz[j][i - (2*int(ring_seperation/step))+4][k]))
				cos_theta[j][k] = np.dot(B_vector, np.array([[0], [size[1]], [0]]))/(np.linalg.norm(B_vector)*np.linalg.norm(np.array([[0], [size[1]], [0]])))
				if(np.isnan(cos_theta[j][k])):
					cos_theta[j][k] = 0

		flux[i] += ((np.sum(Bx[:, i - (2*int(ring_seperation/step))+4, :])**2 + np.sum(By[:, i - (2*int(ring_seperation/step))+4, :])**2 + np.sum(Bz[:, i - (2*int(ring_seperation/step))+4, :])**2)**0.5) * (np.pi/area_element)*ring_radius**2 * abs(np.sum(cos_theta))

plt.plot(z_position, flux)
plt.xlabel("z position (cm)")
plt.ylabel("flux (Wb)")
plt.show()
