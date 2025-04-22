import aar_math_applications as ama

constants = ama.gen_random_poly_func()
#constants = []
point = 5

#matrix = ama.make_arbitrary_matrix()
#matrix = ama.make_particular_arbitrary_matrix(4, 5)
matrix = [[1.0, 2.0, 3.0, -1.0], [2.0, 1.0, 1.0, 1.0], [-4.0, 1.0, 3.0, -5.0], [-9.0, 3.0, 8.0, 12.0]]
#matrix = []

#ama.repr_fraction(3.2/63)
#print()

#print("The square root of 3 is: ", ama.newtons_approx(3))
#print()

#print("The hypotenuse of a right triangle with legs 3 and 5 is: ", ama.pythagorean_theorem(3, 5))
#print()

#print("The square root of 3 is: ", ama.sq_sum([[4, 1], [4, 2]]))
#print()

#print("The generalized pigeonhold principle for (1000, 20) is: ",ama.generalized_pigeonhole_principle(1000, 20))
#print()

#print("The function is: ", constants)
#print("The x we want a particular solution for is: ", point)
#print()

#print(ama.repr_poly_func_particular_sol(constants, point))
#print()

#print("The summation of the sequence from i = 3 to n = 5 is:", ama.gen_summation_notation(constants, 3, 5))
#print()

#print("The summation of the sequence from i = 0 to n = 5 is:\n", ama.repr_gen_summation_notation(constants, 0, 5))
#print()

#print("The derivative of the given function is: ", ama.repr_derivative(constants))
#print()

#print("The antiderivative of the given function is: ", ama.antiderivative(constants))
#print()

#print("The antiderivative of the given function is: ", ama.repr_antiderivative(constants))
#print()

#print(ama.repr_definite_integral(constants, 3, 5))
#print()

#print("The area under the curve of the given function from 3 to 5 is: ", ama.definite_integral(constants, 3, 5))
#print(matrix)

#print("start matrix:\n", ama.repr_matrix(matrix))
#matrix = ama.repr_reduced_row_echelon_form(matrix)
#matrix = ama.clean_matrix(matrix)
#print(ama.repr_matrix(matrix))
#print()

#print(ama.next_collatz(622))

#print(ama.get_vector(matrix, 5))

#matrix = ama.elementary_row_op_scale(matrix, 0, 3)

#print("Matrix after scale row 0 by constant 3\n", ama.repr_matrix(matrix))
#matrix = ama.elementary_row_op_swap(matrix, 0, 1)

#print("Swap rows 0 and 1\n", ama.repr_matrix(matrix))
#matrix = ama.elementary_row_op_elim(matrix, 1, 2, -4)

#print("Use elimination operation: R1 - 4R2\n", ama.repr_matrix(matrix))