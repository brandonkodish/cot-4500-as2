import numpy as np

np.set_printoptions(precision=7, suppress=True)

# Given data points
x = [3.6, 3.8, 3.9]
f = [1.675, 1.436, 1.318]


def neville_method(x_points, y_points, x):
    n = len(x_points)
    Q = [[0.0 for _ in range(n)] for _ in range(n)]

    # Initialize the first column of Q with y-values
    for i in range(n):
        Q[i][0] = y_points[i]

    # Compute the Neville's table
    for i in range(1, n):
        for j in range(n - i):
            xi, xj = x_points[j], x_points[i + j]
            Q[j][i] = ((x - xi) * Q[j + 1][i - 1] - (x - xj) * Q[j][i - 1]) / (xj - xi)

    # The top-right corner of Q contains the interpolated value
    return Q[0][n - 1]


# Interpolate at x = 3.7
x_to_interpolate = 3.7
interpolated_value = neville_method(x, f, x_to_interpolate)
print(" Question 1: ")

print(f"The interpolated value at x = {x_to_interpolate} is {interpolated_value}")

print(" Question 2: \n")


def divided_difference_table(x_points, y_points):
    size = len(x_points)
    matrix = np.zeros((size, size))
    list = []

    for i in range(size):
        matrix[i][0] = y_points[i]

    for i in range(1, size):
        for j in range(1, i + 1):
            matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (x_points[i] - x_points[i - j])

            if i == j:
                list.append(matrix[i][j])

    print(list)
    return matrix


def get_approximate_result(matrix, x_points, value, start):
    reoccuring_x_span = 1
    reoccuring_px_result = start

    for index in range(1, len(matrix)):
        polynomial_coefficient = matrix[index][index]

        reoccuring_x_span *= (value - x_points[index - 1])

        mlt_operation = polynomial_coefficient * reoccuring_x_span

        reoccuring_px_result += mlt_operation

    print(reoccuring_px_result)

x2 = [7.2, 7.4, 7.5, 7.6]
y2 = [23.5492, 25.3913, 26.8224, 27.4589]
divided_table = divided_difference_table(x2, y2)
print(" \n 9.210500000000001, 17.00166666666675, -141.82916666666722 ")

print(" Question 3: \n")

get_approximate_result(divided_table, x2, 7.3, y2[0])
print()

print(" Question 4: \n")

# Given data
x = [3.6, 3.8, 3.9]
f_x = [1.675, 1.436, 1.318]
fp_x = [-1.195, -1.188, -1.182]

n = len(x)
m = 2 * n

# Step 1: Initialize Z (nodes) and F (function values)
Z = []
F = []
df = []

for i in range(n):
    Z.append(x[i])
    Z.append(x[i])
    F.append(f_x[i])
    F.append(f_x[i])
    df.append(fp_x[i])
    df.append(fp_x[i])

# Step 2: Initialize the divided difference table
D = np.zeros((m, m))
for i in range(m):
    D[i][0] = F[i]

# Step 3: Compute first-order divided differences
for i in range(1, m):
    if Z[i] == Z[i - 1]:
        D[i][1] = df[i]
    else:
        D[i][1] = (D[i][0] - D[i - 1][0]) / (Z[i] - Z[i - 1])

# Step 4: Compute higher-order divided differences
for j in range(2, m):
    for i in range(j, m):
        if Z[i] == Z[i - j]:
            # Avoid division by zero for repeated nodes
            D[i][j] = D[i - 1][j - 1] / j
        else:
            D[i][j] = (D[i][j - 1] - D[i - 1][j - 1]) / (Z[i] - Z[i - j])

# Step 5: Print the divided difference table
print("Divided Difference Table:")
for i in range(m):
    for j in range(i + 1):
        print('{0: .6f}'.format(D[i][j]), end='\t')
    print('')

print(" Question 5: \n")

def cubic_spline_interpolation(x_points, y_points):
    size = len(x_points)
    matrix = np.zeros((size, size))
    matrix[0][0] = matrix[size - 1][size - 1] = 1

    for i in range(1, size - 1):
        index = i - 1
        for j in range(index, size, 2):
            matrix[i][j] = x_points[index + 1] - x_points[index]
            index += 1

    for i in range(1, size - 1):
        matrix[i][i] = 2 * ((x_points[i + 1] - x_points[i]) + (x_points[i] - x_points[i - 1]))

    print(np.matrix(matrix), "\n")

    spline_condition = np.zeros((size))

    for i in range (1, size - 1):
        first_term = (3 / (x_points[i + 1] - x_points[i])) * (y_points[i + 1] - y_points[i])
        second_term = (3 / (x_points[i] - x_points[i - 1])) * (y_points[i] - y_points[i - 1])
        spline_condition[i] = first_term - second_term

    print(np.array(spline_condition), "\n")
    print(np.array(np.linalg.solve(matrix, spline_condition)))

x_points, y_points = [2, 5, 8, 10], [3, 5, 7, 9]
cubic_spline_interpolation(x_points, y_points)
