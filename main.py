import math

K = 1  # диэлектрическая постоянная
PHI = 4.5  # Локальный выход электронов
EF = 5.71  # Уровень Ферми
N = 1000


def calculate_j(x, y, z, u):
    z = calculate_z(x, y, z)
    s1 = 3 / (K * PHI)
    s2 = z * (1 - (23 / (((3 * PHI * K * z) + 10) - (2 * u * K * z)))) + s1
    j = 1620 * u * EF * math.exp(-1.025 * z * math.sqrt(calculate_phi_of_z(z, s1, s2, u)))
    return j


def calculate_z(x, y, z):
    return math.sqrt(math.pow(z, 2) + math.pow(x, 2) + math.pow(y, 2))


def calculate_phi_of_z(z, s1, s2, u):
    a = (u * (s1 + s2)) / (2 * z)
    b = 2.86 / (K * (s2 - s1))
    c = math.log((s2 * (z - s1)) / (s1 * (z - s2)))
    result = PHI - a - (b * c)
    return result


def monte_carlo_method():
    pass


if __name__ == '__main__':
    print(calculate_j(2, 2, 2, 0.01))