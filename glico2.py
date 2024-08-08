import math

class Player:
    _tau = 0.5  # System constant that controls the change in volatility
    epsilon = 0.000001

    __q = math.log(10) / 400
    def __init__(self, rating=1500, RD=350, vol=0.06):
        self.rating = rating
        self.RD = RD
        self.vol = vol

    def get_glicko2_rating(self):
        # Convert the rating and RD to Glicko-2 scale
        mu = (self.rating - 1500) * self.__q
        phi = self.RD * self.__q
        return mu, phi
    def get_rating(self):
        return (self.rating / self.__q) + 1500

    def set_rating_deviation(self, rd):
        self.RD = rd / self.__q


    def update_rating(self, opponents, results):
        mu, phi = self.get_glicko2_rating()
        v = self._compute_v(mu, opponents)
        delta = self._compute_delta(mu, opponents, results, v)
        sigma_prime = self._update_volatility(phi, v, delta)
        phi_star = math.sqrt(phi**2 + sigma_prime**2)
        phi_prime = 1 / math.sqrt(1 / phi_star**2 + 1 / v)
        mu_prime = mu + phi_prime**2 * sum(
            self._g(phij) * (sj - self._E(mu, muj, phij))
            for (muj, phij), sj in zip(opponents, results)
        )
        self.RD = phi_prime / self.__q
        self.rating = mu_prime / self.__q + 1500
        self.vol = sigma_prime

    def _compute_v(self, mu, opponents):
        return 1 / sum(
            (self._g(phij) ** 2) * e * (1 - e)
            for muj, phij in opponents
            for e in [self._E(mu, muj, phij)]
        )

    def _compute_delta(self, mu, opponents, results, v):
        return v * sum(
            self._g(phij) * (sj - self._E( mu, muj, phij))
            for (muj, phij), sj in zip(opponents, results)
        )

    def _update_volatility(self, phi, v, delta):
        a = math.log(self.vol ** 2)
        epsilon = 0.000001
        A = a
        if delta ** 2 > phi ** 2 + v:
            B = math.log(delta ** 2 - phi ** 2 - v)
        else:
            k = 1
            while True:
                B = a - k * Player._tau
                if self._f(B, delta, phi, v, a) >= 0:
                    break
                k += 1
        fA = self._f(A, delta, phi, v, a)
        fB = self._f(B, delta, phi, v, a)



        while abs(B - A) > epsilon:
            C = A + (A - B) * fA / (fB - fA)
            fC = self._f(C, delta, phi, v, a)
            if fC * fB < 0:
                A = B
                fA = fB
            else:
                fA /= 2
            B = C
            fB = fC

        return math.exp(A / 2)

    def _f(self, x, delta, phi, v, a):
        exp_x = math.exp(x)
        return exp_x * (delta ** 2 - phi ** 2 - v - exp_x) / 2 / (phi ** 2 + v + exp_x) ** 2 - (x - a) / (Player._tau ** 2)

    def _g(self,phi):
        return 1 / math.sqrt(1 + 3 * phi ** 2 / math.pi ** 2)

    def _E(self, mu, muj, phij):
        return 1 / (1 + math.exp(self._g(phij) * (mu - muj)))


# Example usage
if __name__ == "__main__":
    player = Player(rating=1500, RD=200, vol=0.06)
    opponents = [
        ((1400 - 1500) / 173.7178, 30 / 173.7178),
        ((1550 - 1500) / 173.7178, 100 / 173.7178),
        ((1700 - 1500) / 173.7178, 300 / 173.7178)
    ]
    results = [1, 0, 0]

    player.update_rating(opponents, results)
    print(f"New rating: {player.rating:.10f}")
    print(f"New RD: {player.RD:.10f}")
    print(f"New volatility: {player.vol:.6f}")
