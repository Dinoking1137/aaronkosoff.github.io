import random
#import matplotlib.pyplot as plt
import numpy as np
"""
aar_math_applications.py
==================

This module contains some useful functions for a variety of mathematical applications and for understanding code functions.
This module demonstrates the structure code so that 'help()' calls provide information for each function or class.

USAGE:
======
1) Import the module:
	import aar_math_applications as ama

2) Call functions:
	ama.

3) For detailed help on any function:
	help(ama.)

CONTENTS:
=========
- Constants:
	pi
	num_euler

- Functions:
	Ordering:
		reverse
	
	Algebra:
		newtons_approx

	Geometry:
		2D:
			circumference
			circ_area
			concave_polygonal_perimeter
			convex_polygonal_perimeter
			concave_polygonal_area
			convex_polygonal_area
		
		3D:
			surface_area_sphere
			volume_sphere
			surface_area_cylinder
			volume_cylinder
			surface_area_pyramid_tri
			volume_pyramid_tri
			surface_area_pyramid_sq
			volume_pyramid_sq
			surface_area_pyramid_poly
			volume_pyramid_poly
			surface_area_cone
			volume_cone

	Trigonometry:

		Trig Functions:
			sin
			cos
			tan
			csc
			sec
			cot

		Hyperbolic Trig Functions:
			sinh
			cosh
			tanh
			csch
			sech
			coth

		Arc Trig Functions:
			arc_sin
			arc_cos
			arc_tan
			arc_csc
			arc_sec
			arc_cot

		Arc Hyperbolic Trig Functions:
			arc_sinh
			arc_cosh
			arc_tanh
			arc_csch
			arc_sech
			arc_coth

	Statistics:

	Discrete Mathematics:
		poly_func_particular_sol
		factorial
		permutation
		combination
		generalized_pigeonhole_principle
		summation_notation
		product_notation

	Continuous Mathematics:
		Calculus:
			derivative
			antiderivative
			definite_integral

	Number Theory:
		next_collatz

	Other:
		count_digraph

- Subroutines:
	Pre-Calc:
		repr_elipse
		repr_circle
		repr_hyperbola

	Discrete Mathematics:
		repr_summation

	Continuous Mathematics:
		repr_derivative
		repr_antiderivative
		repr_indefinite_integral

	Definitions:
		# This section will describe what the respective functions actually represent in a mathematical context.
		# The purpose of this section is to educate users on the applications of the given functions past the provided help sections.

Author: Aaron Kosoff
Date: 2025-04-10
"""

#--------------------------------------------------------
# CONSTANTS

pi = 3.14159265359
num_euler = 2.71828182846

#--------------------------------------------------------
# FUNCTIONS

def reverse(arr):
	"""
	Input an array
	The outputted array, string, or manipulable type will be of reverse order
	"""
	return arr[::-1]

def gen_random_poly_func():
	""""""
	arr_len = random.randint(1, 10)
	arr_poly = []

	for i in range(arr_len):
		arr_poly += [random.randint(1, 10)]
		i+=1

	return arr_poly

def spec_gen_random_poly_func(l):
	""""""
	arr_poly = []

	for i in range(l):
		arr_poly += [random.randint(1, 10)]
		i+=1

	return arr_poly

def make_arbitrary_matrix():
	"""
	Creates an arbitrary k*n matrix where both k and n are randomly generated from the range [1, 10].
	After generating the rows and columns, the function then sets random values in each space from the range [-10, 10].
	"""
	row_len = random.randint(1, 10)
	col_len = random.randint(1, 10)
	matrix = []

	for r in range(row_len):
		matrix.insert(r, [])
		for c in range(col_len):
			matrix[r].append(random.randint(-10, 10))

	return matrix

def make_particular_arbitrary_matrix(row_len, col_len):
	"""
	Creates an arbitrary k*n matrix where both k and n are set in the function call.
	The row length and column length both have to be greater than or equal to 0 for this function to work properly.
	After generating the rows and columns, the function then sets random values in each space from the range [-10, 10].
	When making, adding an extra column will give the appended B vector for the linear system.
	"""

	matrix = []

	for c in range(col_len):
		matrix.insert(c, [])
		for r in range(row_len):
			matrix[c].append(random.randint(-10, 10))

	return matrix

def newtons_approx(sq_num):
	"""
	This function takes in a square as the argument and will find its square under a very small room for error
	In its current state it works best with numbers >= 1
	"""
	if(sq_num < 0):
		return 0

	sqr = sq_num

	while(True):
		if(abs(((sqr**2) - sq_num)) < 0.00001):
			break
		sqr = (sqr + (sq_num / sqr)) / 2.0;

	return sqr

def sq_sum(sq_arr):
	"""
	This function takes in the 2d array of squares, where each element of the array is assumed to be an embedded array.
	The function will take the first value of the sq_array and take it as the square that is being multiplied by the given constants.
	The consecutive values will be taken as constants that are multiplied with the first value, and then added to the output.
	"""
	#doesnt work
	out_sq = 0.0

	for i in range(len(sq_arr)):
		if (len(sq_arr[i]) > 1):
			for j in range(len(sq_arr[i-1])):
				out_sq += sq_arr[i, j+1] * newtons_approx(sq_arr[i, 0])
			break
		out_sq += newtons_approx(sq_arr[i])

	return out_sq

def degrees_to_radians(d):
	"""
	Takes in the specified degree value to convert to radians.
	"""
	return d * (pi/180.0)

def radians_to_degrees(r):
	"""
	Takes in a radian value to convert to degrees.
	The value given should be a fraction of the given pi value from the library.
	The value inputed should be of the form (c * ama.pi)/n where c and n are particular constants.
	"""
	return r * (180.0/pi)

def circumference(r):
	"""
	Returns the circumference of a circle for the given radius.
	"""
	return 2.0 * pi * r

def circ_area(r):
	"""
	Returns the radius of a circle for the given radius.
	"""
	return pi * (r**2)

def pythagorean_theorem(a, b):
	"""
	Returns c for the equation a^2 + b^2 = c^2
	"""
	c = (a**2) + (b**2)
	return newtons_approx(c)

def concave_polygonal_perimeter(r, sides):
	"""
	CAUTION: THIS FUNCTION ASSUMES THE GIVEN POLYGON IS BOTH EQUIANGULAR AND EQUILATERAL
	Returns the perimeter for the given "radius".
	The "radius" for this formula should be radius of the inscribed circle.
	"""
	# = 6(2r/sqrt(3))
	#work with the given function to find a universal one
	#The function presented is for finding the perimiter of a hexagon

def convex_polygonal_perimeter(r_outer, r_inner, sides):
	"""
	CAUTION: THIS FUNCTION ASSUMES THE GIVEN POLYGON IS BOTH EQUIANGULAR AND EQUILATERAL
	Returns the perimeter for the given radii.
	The variable r_outer is representative of the line that extends from the center of the polygon connecting to one of the outer vertices.
	The variable r_inner is representative of the line that extends from the center of the polygon connecting to one of the inner vertices which wouldn't exist if the function were concave.
	"""

def concave_polygonal_area(r, sides):
	"""
	CAUTION: THIS FUNCTION ASSUMES THE GIVEN POLYGON IS BOTH EQUIANGULAR AND EQUILATERAL
	Returns the perimeter for the given "radius".
	The "radius" for this formula should be radius of the inscribed circle.
	"""
	# A = (sides/2) * (side length) * (r)

def convex_polygonal_area(r_outer, r_inner, sides):
	"""
	CAUTION: THIS FUNCTION ASSUMES THE GIVEN POLYGON IS BOTH EQUIANGULAR AND EQUILATERAL
	Returns the perimeter for the given radii.
	The variable r_outer is representative of the line that extends from the center of the polygon connecting to one of the outer vertices.
	The variable r_inner is representative of the line that extends from the center of the polygon connecting to one of the inner vertices which wouldn't exist if the function were concave.
	"""
	# concave_polygonal_area - ((sides/2) * ((r_outer - r_inner) * (Distance between 2 vertices /2))

def poly_func_particular_sol(c, x_position):
	"""
	Define an array in the order of the constants that are attached to variables in the function.
	The function should be ordered by power
	The set of constants inputted are assumed to be appended to these variables in this set order
	"""
	output = 0.0

	for i in range(len(c)):
		output += ((x_position)**i) * c[i]

	return output

def factorial(n):
	"""
	Info
	"""
	if(n<=1):
		return 1
	return n*factorial(n-1)

def permutation(n,r):
	"""
	Info
	"""
	return factorial(n)/factorial(n-r)

def combination(n, r):
	"""
	Info
	"""
	return permutation(n, r)/factorial(r)

def generalized_pigeonhole_principle(k, b):
	"""
	The pigeonhole principle describes the minimum value where a certain condition HAS to happen
	K is representative of the size of the size of the domain
	B is representative of the goal
	"""
	return (int)((k * b) - k + 1)

def gen_summation_notation(c, i, n, cycle_output = 0):
	"""
	Returns the value generated from a summation function.
	Takes c as the array of summation terms.
	Takes n as the upper bound of the summation.
	Takes i as the lower bound of the summation and the increment value.
	Cycle output is just what is returned when i = n.
	"""

	for x in range(len(c)):
		cycle_output += ((i)**(x)) * c[x]

	if (i<n):
		return gen_summation_notation(c, i+1, n, cycle_output)

	return cycle_output

def antiderivative(c):
	"""
	Returns the antiderivative created for a polynomial function with the coefficients defined in an array cast as the argumentative call.
	"""
	terms = f""
	c = reverse(c)
	c_return = ["c"]

	for i in range(len(c)):
		c_return = [c[i]/(i+1)] + c_return

	return c_return

def definite_integral(c, a, b):
	"""
	This function takes the given polynomial function and the ranges for the integral and returns a float value representative of the area under the curve.
	"""

	c = reverse(c)
	out_a = 0.0
	out_b = 0.0

	for i in range(len(c)):
		out_a += (c[i]/(i+1)) * ((a)**(i+1))
		out_b += (c[i]/(i+1)) * ((b)**(i+1))

	return (out_b - out_a)

def next_collatz(n):
	"""
	Obtains the next collatz number.
	If the argument value is even, then the function will return n/2
		Checks parity with the conditional statement if(n%2 == 0)
	If the argument value is odd, then the function will return 3n + 1
		Checks parity with the conditional statement if(n%2 == 1)
	If neither conditionals are checked, the function returns "ERROR"
	"""
	if(n%2 == 0):
		return (int)(n/2)
	if(n%2 == 1):
		return (3*n) + 1
	return f"ERROR"

#--------------------------------------------------------
# ROUTINES

def repr_fraction(num, check=1):
	"""
	Idea:
	Checks incrementally to see what the denominator must've been of a number, and therefore the numerator by checking the error between the value and the int version of the value for 0
	"""
	check_num = num * check

	if(check_num > (int)(check_num)):
		#print(check_num, " > ", (int)(check_num))
		repr_fraction(num, check+1)
	else:
		#print(check_num, " <= ", (int)(check_num))
		print((int)(check_num), "/", (check))

def repr_fraction_addition(num_1, num_2):
	""""""

def repr_degrees_to_radians(d):
	"""
	Takes in the specified degree value to convert to radians.
	"""
	return f"{d/180.0}pi"

def repr_func(c):
	"""
	Prints the given polynomial function with the coefficients defined in an array cast as the argumentative call as a format string.
	"""
	terms = f""

	for i in range(len(c)):
		if(i > 0):
			terms = f"{c[i]}x^{i} + " + terms
		else:
			terms = f"{c[i]}" + terms

	return terms

def points_of_intersect(c1, c2):
	"""
	Finds the points of intersection for functions f1 and f2.
	The arguments c1 and c2 represent the set of coefficients for the function.
	The coefficients will automatically be set to the respective x terms.
	"""

def repr_poly_func_particular_sol(c, x_position):
	"""
	This function takes in the given polynomial function and the given x value and finds the particular solution at that said x value.
	The function returns the process of finding the particular solution, and finally represents the answer within the last step.
	The function returns these steps through a formating string.
	"""
	c = reverse(c)
	equation = f"For the function {repr_func(c)}, \nthe particular solution at x = {x_position} is found through the following steps: \n\n= "
	terms = f""
	simp_terms = f""

	for i in range(len(c)):

		if(i > 0):
			terms = f"{c[i]}({x_position})^{i} + " + terms
			simp_terms = f"{((x_position)**i)*c[i]} + " + simp_terms

		else:
			terms = f"{c[i]}({x_position})^{i} \n= " + terms
			simp_terms = f"{((x_position)**i)*c[i]} \n= " + simp_terms

	terms += simp_terms
	terms += f"{poly_func_particular_sol(c, x_position)}"
	equation += terms
	return equation

def repr_table(c, start, end, step = 2):
	"""
	This function returns the table of the given equation c.
		The c argument is to be defined
	The start represents the initial value the table uses.
	The end represents the final value that the table will reach.

	The step value is defaulted to 1, however it can be redifined to show how many elements in the table there will be.
	The function will automatically convert this value to int.
	The step is the amount of elements, not the size x increments within the table.
	"""

def ellipse():
	""""""

def repr_hyperbola():
	""""""

def polar_function():
	"""
	Defines a function in polar coordinates.
	Potential visual element to come in the future.
	"""

def repr_gen_summation_notation(c, i, n):
	"""
	Returns the format string detailing the summation from i to n, detailing every step present.
	"""
	#Work on this one
	terms = f""
	c = reverse(c)

	for y in range(n-(i-1)):
		out = 0.0
		terms_set = f""

		for x in range(len(c)):
			if(x > 0):
				out += c[x]*((y+i)**(x))
				terms_set = f"{c[x]*((y+i)**(x))} + " + terms_set
			else:
				out += c[x]*((y+i)**(x))
				terms_set = f"{c[x]*((y+i)**(x))} = " + terms_set

		terms_set += f"{out}\n"
		terms_set = f"i = {i+y}: " + terms_set
		terms = terms + terms_set

	terms += f"= {gen_summation_notation(reverse(c), i, n)}"

	return terms

def repr_inst_roc():
	"""
	Returns a string detailing the fundamental definitions behind the instantaneous rate of change.
	"""
	return 0

def repr_concept_poly_derivative():
	""""""
	return 0

def repr_derivative(c):
	"""
	Creates the derivative for a polynomial function with the coefficients defined in an array cast as the argumentative call.
	"""
	terms = f""
	c = reverse(c)

	for i in range(len(c)):
		if(i > 1):
			terms = f"{c[i]*i}x^{i-1} + " + terms
		elif(i == 1):
			terms = f"{c[i]*i}" + terms
		else:
			terms += f" + 0"

	return terms

def repr_antiderivative(c):
	"""
	Creates the antiderivative for a polynomial function with the coefficients defined in an array cast as the argumentative call.
	"""
	terms = f""
	c = reverse(c)

	for i in range(len(c)):
		if (i > 0):
			terms = f"{c[i]}/{i+1}x^{i+1} + " + terms
		else:
			terms = f"{c[i]}x + " + terms

	terms += f"C"

	return terms

def repr_definite_integral(c, a, b):
	"""
	This function takes the given polynomial function, and the ranges of the integral to be taken.
	The function returns the process of finding the area under the curve through a variety of steps and is returned within a format string.
	"""

	equation = f"For the function {repr_func(c)}, \nthe definite integral from {a} to {b} is found through the following steps: \n\n= "

	equation += repr_antiderivative(c)
	equation += "\n= "

	c = reverse(c)
	out_a = 0.0
	out_b = 0.0

	repr_out_a = f""
	repr_out_b = f""

	repr_out_a_simp = f""
	repr_out_b_simp = f""

	for i in range(len(c)):
		if(i > 0):
			repr_out_a = f"{c[i]}/{(i+1)}({a})^{(i+1)} + " + repr_out_a
			repr_out_a_simp = f"{(c[i]/(i+1)) * ((a)**(i+1))} + " + repr_out_a_simp

			repr_out_b = f"{c[i]}/{(i+1)}({b})^{(i+1)} + " + repr_out_b
			repr_out_b_simp = f"{(c[i]/(i+1)) * ((b)**(i+1))} + "+ repr_out_b_simp

		else:
			repr_out_a = f"{c[i]}/{(i+1)}({a})^{(i+1)}) \n= " + repr_out_a
			repr_out_a_simp = f"{(c[i]/(i+1)) * ((a)**(i+1))}) \n= " + repr_out_a_simp

			repr_out_b = f"{c[i]}/{(i+1)}({b})^{(i+1)})" + repr_out_b
			repr_out_b_simp = f"{(c[i]/(i+1)) * ((b)**(i+1))}" + repr_out_b_simp

	equation += f"(" + repr_out_b + f") - ("
	equation += repr_out_a

	equation += f"(" + repr_out_b_simp + f") - ("
	equation += repr_out_a_simp

	equation += f"{definite_integral(c, a, b)}"

	return equation

def repr_area_between_curves(c1, c2, a, b):
	"""
	Will find the area between two specified curves.
	Will find points of intersection and determine places where either function is greater than the other.
	Follows the formula fn_int(f1, a, b) - fn_int(f2, a, b) for all places where f1 is greater than f2.
	"""

def taylor_series(c):
	"""
	Make taylor series of something lmao
	"""
	#not implemented
	return 0

def polynomial_interpolation():
	"""
	Wimmel Code.
	"""
	return 0

def lagrange_interpolation():
	"""
	Wimmel Code.
	"""
	return 0

def fibonacci(n):
	""""""

def check_matrix(m):
	"""
	Checks if the given matrix is valid, I.E. each row should have the same amount of column elements.
	Matrices are valid even if the amount of rows isn't equal to the amount of columns.
	Matrices are valid even if it is empty.
	"""
	check = [0, 0]
	for row in range(len(m)-1):
		print()

def repr_matrix(m):
	"""
	Represents a matrix.
	Vectors should be added as an embedded array.
	Can represent any matrix, inserting an empty array will give a matrix of zero rows and zero columns.
	"""
	matrix_out = f""

	if not m or not m[0]:
		return matrix_out

	num_rows = len(m[0])
	num_cols = len(m)

	for row in range(num_rows):
		for col in range(num_cols):
			matrix_out += f"{m[col][row]:.3f} "
		matrix_out += "\n"
	return matrix_out

def append_vector(m, v):
	"""
	Will append a vector to the end of a matrix.
	Used ideally for applying the vector of the set b values.
	Can even be done for an empty matrix.
	"""

def clean_matrix(m):
	""""""
	num_rows = len(m[0])
	num_cols = len(m)

	for row in range(num_rows):
		for col in range(num_cols):
			m[col][row] = repr_fraction(m[col][row])

def get_vector(m, vect_num):
	"""
	Returns the given vector if present in the matrix.
	Otherwise will return a statement that indicates that the vector is out of bounds.
	"""
	out = f"The vector is: \n"
	if(vect_num <= len(m)):
		for r in range(len(m[vect_num])):
			out += f"{m[vect_num][r]:.3f}\n"
	else:
		out += f"ERROR: OUT OF BOUNDS"
	return out

def elementary_row_op_scale(m, v, c):
	"""
	Follows the formula: m[v] -> c*m[v] 
	Performs the scale elementary function on a matrix.
	The scalar and vector are specified on call.
	v is simply the index value of the vector, not the vector itself.
	Does not alter the values represented within the given linear equation.
	"""
	num_cols = len(m)

	for col in range(num_cols):
		m[col][v] *= c

	return m

def elementary_row_op_swap(m, v, u):
	"""
	Follows the formula: m[v] <=> m[u]
	Performs the swap elementary function on a matrix.
	The vectors to be swapped are specified on call.
	v and u are simply the index values of their vectors, not the vectors themselves.
	Does not alter the values represented within the given linear equation.
	"""
	num_cols = len(m)

	for col in range(num_cols):
		temp = m[col][v]
		m[col][v] = m[col][u]
		m[col][u] = temp

	return m

def elementary_row_op_elim(m, v, u, c):
	"""
	Follows the formula: m[v] -> m[v] + c*m[u]
	Performs the elimination elementary function on a matrix.
	The scalar and vectors are specified on call.
	v and u are simply the index values of their vectors, not the vectors themselves.
	Does not alter the values represented within the given linear equation.
	"""
	num_cols = len(m)

	for col in range(num_cols):
		m[col][v] += (c*m[col][u])

	return m

def gaussian_elimination(m, start_row = 0, start_col = 0, pivot_col = False):
	"""
	Performs all elementary row operations on the given matrix through the gaussian elimination process.
	This version continues until it reaches row echelon form.
	To get the b values, simply add an extra b vector, a new program to append b vectors will be added shortly.
	1. 
		Beginning on the left, find the first column that has a nonzero entry. We refer to this as the jth column
	2. 
		A. If the first row of the jth column has a zero entry, obtain a nonzero entry in the first row by exchanging the first row with any row that has a nonzero entry
		B. Use a scaling operation to make the entry in the first row 1
	3. 
		For each i >= 2, make each entry aij equal to zero by applying the elimination operation: Ri-aijR1+Ri
	4. 
		Ignore the top row. Repeat this process with all the rows
			Stopping at step 4 gives the echelon form
	5. 
		Starting with the leading one in the last nonzero row work upwards:
			For each non-zero entry above the leading one make it zero by applying the elimination operation

	# 1 says check from start_row for length of column till a nonzero is found
	# 2.A says elementary_row_op_swap(m, row, start_row)
	# 2.B says elementary_row_op_scale(m, start_row, 1/m[row][start_row])
	# 3 says for all other rows elementary_row_op_elim(m, v, start_row, 1/m[row][v])
	# 4 implies if statement into recursion, where gaussian_elimination(m, start_row + 1, start_col + 1)
		# This added start_row is done to keep the process going for the entire matrix
	# 5 now finds the pivot columns and applies a similar elimination tactic as #3 for all columns above it
		# No below since the below should already be removed by 3
		# If there are still terms under the pivot columns, then there is an error

	"""

	#IMPORTANT: m[c][r]

	if(m):
		for r in range(start_row, len(m[0])):
			if(m[start_col][r] != 0):

				pivot_col = True

				print(f"R{start_row} <--> R{r}\n")
				m = elementary_row_op_swap(m, start_row, r)
				
				print(f"R{start_row} -> ({1.0/(m[start_col][start_row]):.3f})R{r}\n")
				m = elementary_row_op_scale(m, start_row, 1.0/(m[start_col][start_row]))
				print(repr_matrix(m))
				break

		for r in range(start_row + 1, len(m[0])):

			print(f"R{r} -> R{r} - ({m[start_col][r]:.3f})R{start_row}\n")
			elementary_row_op_elim(m, r, start_row, -m[start_col][r])

		print(repr_matrix(m))

		if(pivot_col and (start_row < len(m[0])-1) and (start_col < len(m)-1)):
			gaussian_elimination(m, start_row+1, start_col+1)

		if(pivot_col is False and (start_row < len(m[0])-1) and (start_col < len(m)-1)):
			gaussian_elimination(m, start_row, start_col+1)

	return m

def repr_reduced_row_echelon_form(m):
	"""
	Performs all elementary row operations on the given matrix through the gaussian elimination process.
	This version continues until it reaches reduced row echelon form.
	"""
	m = gaussian_elimination(m)

def repr_pivot_free(m):
	"""
	This function returns an f-string that details which columns are pivot and which ones are free.
	To run this program, please run the gaussian elimination function as otherwise it will return inaccurate results.

	There may be a safety measure in place for this, however please ensure safe coding practice.
	"""

def dot_product(u, v):
	"""

	"""

def check_matrix_compatibility(mA, mB):
	"""
	For this function, have mA defined by k*n and mB defined by l*k.
	The ammount of rows mB should match the amount of columns mA.
	This means that the multiplied version of the matrix will be a l*n matrix.

	This function does not require that the standards for a compatible matrix be met.
	This function simply returns true if applicable and false otherwise.
	
	RUN THIS FUNCTION PRIOR TO RUNNING THE MATRIX MULTIPLICATION FUNCTION!!!
	"""


def matrix_multiplication(mA, mB):
	"""
	For this function, have mA defined by k*n and mB defined by l*k.
	The ammount of rows mB should match the amount of columns mA.
	This means that the multiplied version of the matrix will be a l*n matrix.

	This function REQUIRES that the standards for a compatible matrix be met.
	IF the functions are not met, unknown errors ensue.
	"""

def collatz_steps_to_convergence(n):
	"""
	Finds how many steps are taken for a number to converge to 1 using the collatz conjecture.
	Be careful not to implement a number thats too large for the system.
	If dealing with a sufficiently large number, the system will experience what is known as a "Recursion Limit".
	This limit usually occurs with values that need to take ~1000 steps to converge, so please be careful!
	"""