import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import lu_factor, lu_solve


def interpolate(points, x):
    result = 0

    for i in range(len(points)):
        xi, yi = points[i]
        base = 1

        for j in range(len(points)):
            if i == j:
                continue

            xj, yj = points[j]
            base *= (float(x) - float(xj)) / (float(xi) - float(xj))

        result += float(yi) * base

    return result


def interpolate_with_lagrange(k_values):
    for filename in os.listdir('./topografia'):
        file_path = './topografia/' + filename
        with open(file_path, 'r') as f:
            topografia = list(csv.reader(f))

        fig, axs = plt.subplots(1, len(k_values), figsize=(6 * len(k_values), 5))

        for i, k in enumerate(k_values):
            interpolation = topografia[1::k]



            interpolated_height = []
            distance = []
            height = []

            for point in topografia[1:]:
                x, y = point

                distance.append(float(x))
                height.append(float(y))
                interpolated_height.append(interpolate(interpolation, float(x)))

            train_distance = []
            train_height = []
            for point in interpolation:
                x, y = point
                train_distance.append(float(x))
                train_height.append(interpolate(interpolation, float(x)))

            axs[i].semilogy(distance, height, 'r.', label='pełne dane')
            axs[i].semilogy(distance, interpolated_height, color='blue', label='funkcja interpolująca')
            axs[i].plot(train_distance, train_height, 'g.', label='dane do interpolacji')
            axs[i].legend()
            axs[i].set_ylabel('Wysokość')
            axs[i].set_xlabel('Odległość')
            axs[i].set_title('Przybliżenie interpolacją Lagrange\'a, ' + str(len(interpolation)) + ' punkty(ów)')
            axs[i].grid()

            axs[i].set_yscale('log')
            min_height = min(filter(lambda x: x > 0, height))
            max_height = max(height)
            axs[i].set_ylim([min_height / 2, max_height * 2])

        fig.suptitle(filename)
        plt.tight_layout()
        plt.show()


def calculate_params(points):
    n = len(points)
    A = np.zeros((4 * (n - 1), 4 * (n - 1)))
    b = np.zeros(4 * (n - 1))

    for i in range(n - 1):
        x, y = points[i]
        A[4 * i + 3, 4 * i + 3] = 1
        b[4 * i + 3] = float(y)

    for i in range(n - 1):
        x1, y1 = points[i + 1]
        x0, y0 = points[i]
        h = float(x1) - float(x0)
        A[4 * i + 2, 4 * i] = h ** 3
        A[4 * i + 2, 4 * i + 1] = h ** 2
        A[4 * i + 2, 4 * i + 2] = h
        A[4 * i + 2, 4 * i + 3] = 1
        b[4 * i + 2] = float(y1)

    for i in range(n - 2):
        x1, y1 = points[i + 1]
        x0, y0 = points[i]
        h = float(x1) - float(x0)
        A[4 * i, 4 * i] = 3 * (h ** 2)
        A[4 * i, 4 * i + 1] = 2 * h
        A[4 * i, 4 * i + 2] = 1
        A[4 * i, 4 * (i + 1) + 2] = -1

    for i in range(n - 2):
        x1, y1 = points[i + 1]
        x0, y0 = points[i]
        h = float(x1) - float(x0)
        A[4 * (i + 1) + 1, 4 * i] = 6 * h
        A[4 * (i + 1) + 1, 4 * i + 1] = 2
        A[4 * (i + 1) + 1, 4 * (i + 1) + 1] = -2

    A[1, 1] = 2
    x1, y1 = points[-1]
    x0, y0 = points[-2]
    h = float(x1) - float(x0)
    A[-4, 1] = 2
    A[-4, -4] = 6 * h

    return lu_solve(lu_factor(A), b)


def interpolation_function(points, x, params):
    n = len(points)
    param_array = np.reshape(params, (n - 1, 4))

    for i in range(1, len(points)):
        xi, yi = points[i - 1]
        xj, yj = points[i]
        if float(xi) <= float(x) <= float(xj):
            a, b, c, d = param_array[i - 1]
            h = float(x) - float(xi)
            print(i, x, xi, h)
            return a * (h ** 3) + b * (h ** 2) + c * h + d

    return -123


def interpolate_with_spline(k_values):
    for file in os.listdir('./topografia'):
        with open('./topografia/' + file, 'r') as f:
            topografia = list(csv.reader(f))
            topografia = topografia[1:]

            fig, axs = plt.subplots(1, len(k_values), figsize=(6 * len(k_values), 5))

            for i, k in enumerate(k_values):
                shift = (-1) * (len(topografia) % k)
                if shift != 0:
                    interpolation_data = topografia[:shift:k]
                else:
                    interpolation_data = topografia[::k]



                distance = []
                height = []
                interpolated_height = []

                for point in topografia:
                    x, y = point
                    distance.append(float(x))
                    height.append(float(y))
                    params = calculate_params(interpolation_data)
                    interpolated_height.append(interpolation_function(interpolation_data, x, params))

                train_distance = []
                train_height = []

                for point in interpolation_data:

                    x, y = point

                    train_distance.append(float(x))
                    train_height.append(float(y))


                shift = -1 * interpolated_height.count(-123)

                axs[i].plot(distance, height, 'r.', label='pełne dane')
                axs[i].plot(distance[:shift], interpolated_height[:shift], color='blue', label='funkcja interpolująca')
                axs[i].plot(train_distance, train_height, 'g.', label='dane do interpolacji')
                axs[i].legend()
                axs[i].set_ylabel('Wysokość')
                axs[i].set_xlabel('Odległość')
                axs[i].set_title(f'Przybliżenie interpolacją Splajnami, {len(interpolation_data)} punkty(ów)')
                axs[i].grid()

            fig.suptitle(file)
            plt.tight_layout()
            plt.show()


if __name__ == '__main__':
    k_values = [10, 25, 40, 50]  # Lista różnych wartości k
    interpolate_with_lagrange(k_values)
    interpolate_with_spline(k_values)