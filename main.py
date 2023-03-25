import math

K = 1  # диэлектрическая постоянная
PHI = 4.5  # Локальный выход электронов
EF = 5.71  # Уровень Ферми
N = 1000
MAX_RND = 32767
SEED = 65539


class RndMultiCmpGenerator:
    def __init__(self, x):
        self.x = x
        self.m = 32767  # 2147483647
        self.a = 65539

    def next(self):
        self.x = (self.a * self.x) % self.m
        return self.x


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


def monte_carlo_method(xl, xh, yl, yh, z, u, random):
    a = [xl, yl]
    b = [xh, yh]
    integral = 0.0
    v = (xh - xl)*(yh - yl)
    x_buff = [0.0, 0.0]
    flag = True
    for i in range(N):
        for j in range(2):
            x_buff[j] = a[j] + (b[j] - a[j]) * random[i][j]
        flag = True
        for j in range(2):
            if (x_buff[j] < a[j]) and (x_buff[j] > b[j]):
                flag = False
        if flag:
            integral = calculate_j(x_buff[0], x_buff[1], z, u)
    integral = (integral * v) / N
    return integral


if __name__ == '__main__':
    random_generator = RndMultiCmpGenerator(SEED)
    random_values = [[random_generator.next() / MAX_RND for x in range(2)] for y in range(N)]