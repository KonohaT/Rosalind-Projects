import math
import matplotlib.pyplot as plt
import numpy as np


def binomial_plot(x: int, prob: float):
    def binomial_distribution(n_chances, p_probability, r_occurrences):
        binomial_coefficient = math.factorial(n_chances)/(math.factorial(r_occurrences)*math.factorial(n_chances - r_occurrences))
        probability = binomial_coefficient * p_probability**r_occurrences * (1-p_probability)**(n_chances - r_occurrences)
        return probability
    
    y_vals = [binomial_distribution(x, prob, i) for i in range(x)]
    x_vals = [i for i in range(x)]
    _, ax = plt.subplots()
    ax.plot(x_vals, y_vals)
    plt.show()

def inverse_binomial_plot(n_chances, r_occurrences): #given a specific number of occurrences (e.g. 30% of 64 offsprings being AaBb), what are the most likely odds of AaBb appearing?
    def binomial_distribution(n_chances, p_probability, r_occurrences):
        binomial_coefficient = math.factorial(n_chances)/(math.factorial(r_occurrences)*math.factorial(n_chances - r_occurrences))
        probability = binomial_coefficient * p_probability**r_occurrences * (1-p_probability)**(n_chances - r_occurrences)
        return probability
    
    y_vals = [binomial_distribution(n_chances=n_chances, p_probability=prob, r_occurrences=r_occurrences) for prob in np.arange(0, 1, 0.01)]
    x_vals = [prob for prob in np.arange(0, 1, 0.01)]
    _, ax = plt.subplots()
    ax.plot(x_vals, y_vals)
    plt.show()

#The odds that a lame antelope (Steve) is caught, given the fraction of overall deer that are caught (caught_fraction)and the fraction of caught deer that are lame (lame_in_caught fraction) or healthy.
#probability of caught given lame; caught = librarian, lame = meek, healthy = not meek
def bayes_antelope(caught_frac: float, lame_caught_frac: float, lame_uncaught_frac: float): 
    uncaught_frac = 1 - caught_frac
    odds_caught = caught_frac * lame_caught_frac / (caught_frac * lame_caught_frac + uncaught_frac * lame_uncaught_frac)
    print(odds_caught)

#bayes_antelope(caught_frac = 1/21, lame_caught_frac=0.4, lame_uncaught_frac=0.1)

#the odds that Steve is a librarian given that he's meek and that 40% of librarians are meek and 10% of non-librarians are meek, 
#and there is 1 librarian for every 20 non-librarians
def bayes_librarian(meek_librarian_frac = 0.4, meek_nonlib_frac = 0.1, librarian_pop_frac = 1/21): 
    """Implements Bayes Theorum from the 3Blue1Brown example"""
    steve_lib_odds = librarian_pop_frac * meek_librarian_frac / ((librarian_pop_frac * meek_librarian_frac) + (1 - librarian_pop_frac) * meek_nonlib_frac)
    print(steve_lib_odds)

bayes_librarian()
