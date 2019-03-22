# The system_state.py file generates all the steady state probabilities and transition rates for our various systems.
# This is made a separate file, as it is referenced in both the Markov processes and the Monte Carlo simulations.
import numpy as np


# Creates a class for each of the systems we're concerned with.
class System(object):
    def __init__(self, name):
        # Gives a name to the system (ie: G1/G2, G3, etc.)
        self.name = name
        # Stores the states as a dictionary of objects, where each key is the name of the state, and the value
        # of the key is the object.
        self.states = dict()


# Creates a class for each individual state we're concerned with.
class State(object):
    def __init__(self, name):
        # Gives a name to the state (ie: 1, 2, 3,..)
        self.name = name
        # Defines a dictionary for transition rates, where each key is the state being transitioned to,
        # and the value of the key is the rate.
        self.transition = dict()
        # Gives a steady state probability for each state.
        self.ss_prob = 0


# Gets the steady state probabilities and transition rates for our final generator states.
def gen_info():
    # All transition rates will be given in time units of days.
    lambda_g = 0.1
    mu_g = 3
    # Looking at each individual generator first.
    p_gu = mu_g / (mu_g + lambda_g)
    p_gd = lambda_g / (mu_g + lambda_g)
    # Since G3 can be considered a 2 state system, we construct a system for it.
    two_state_gen = System('G3 System')
    for i in range(1, 3):
        temp_state = State(i)
        if i == 1:
            temp_state.transition[2] = lambda_g
            temp_state.ss_prob = p_gu
        else:
            temp_state.transition[1] = mu_g
            temp_state.ss_prob = p_gd
        two_state_gen.states[i] = temp_state
    # Adds in state objects to the system class, along with transition rates and probabilities.
    # This is pretty much hard coding values, but we try to keep this format for more complex systems.
    # Now we look at a four state system, and begin by initializing a system class.
    four_state_gen = System('4 State Generation')
    for j in range(1, 5):
        # Creates a temporary state object.
        temp_state = State(j)
        # Checks for the Up/Up state.
        if j == 1:
            temp_state.transition[2] = lambda_g
            temp_state.transition[3] = lambda_g
            temp_state.ss_prob = np.square(p_gu)
        # Checks for the Down/Down state.
        elif j == 4:
            temp_state.transition[2] = mu_g
            temp_state.transition[3] = mu_g
            temp_state.ss_prob = np.square(p_gd)
        # Checks for the Up/Down or Down/Up states, which are equivalent.
        else:
            temp_state.transition[1] = mu_g
            temp_state.transition[4] = lambda_g
            temp_state.ss_prob = p_gu * p_gd
        four_state_gen.states[j] = temp_state
    # We also reduce this 4 state system to a 3 state one, where the states are sorted by their amount of generation.
    # Calculating the steady state probabilities of each state, and adding them to the dictionary.
    # We'll go ahead and call this the G1/G2 system, since we'll be using it.
    three_state_gen = System('G1/G2 System')
    for k in range(1, 4):
        temp_state = State(k)
        if k == 1:
            temp_state.transition[2] = 2 * lambda_g
            temp_state.ss_prob = four_state_gen.states[1].ss_prob
        elif k == 2:
            temp_state.transition[3] = lambda_g
            temp_state.transition[1] = mu_g
            temp_state.ss_prob = four_state_gen.states[2].ss_prob + four_state_gen.states[3].ss_prob
        else:
            temp_state.transition[2] = 2 * mu_g
            temp_state.ss_prob = four_state_gen.states[4].ss_prob
        three_state_gen.states[k] = temp_state
    # Returns the three state system.
    return two_state_gen, three_state_gen


# Gets the steady state probabilities and transition rates for our final transmission states.
def transmission_info():
    # All transition rates will be given in time units of days.
    lambda_tn = 10 / 365
    lambda_ta = 100 / 365
    mu_t = 3
    n = 0.12
    s = 1.2
    # Calculating the steady state probability of a transmission line being up or down.
    p_tun = mu_t / (mu_t + lambda_tn)
    p_tua = mu_t / (mu_t + lambda_ta)
    p_tdn = lambda_tn / (mu_t + lambda_tn)
    p_tda = lambda_ta / (mu_t + lambda_ta)
    # Creates a four state system for normal weather.
    four_state_n = System('4 State Transmission, Normal')
    for i in range(1, 5):
        temp_state = State(i)
        if i == 1:
            temp_state.transition[2] = lambda_tn
            temp_state.transition[3] = lambda_tn
            temp_state.ss_prob = np.square(p_tun)
        elif i == 4:
            temp_state.transition[2] = mu_t
            temp_state.transition[3] = mu_t
            temp_state.ss_prob = np.square(p_tdn)
        else:
            temp_state.transition[1] = mu_t
            temp_state.transition[4] = lambda_tn
            temp_state.ss_prob = p_tun * p_tdn
        four_state_n.states[i] = temp_state
    # Creates a three state system for normal weather.
    three_state_n = System('3 State Transmission, Normal')
    for j in range(1, 4):
        temp_state = State(j)
        if j == 1:
            temp_state.transition[2] = 2 * lambda_tn
            temp_state.ss_prob = four_state_n.states[1].ss_prob
        elif j == 2:
            temp_state.transition[3] = lambda_tn
            temp_state.transition[1] = mu_t
            temp_state.ss_prob = four_state_n.states[2].ss_prob + four_state_n.states[3].ss_prob
        else:
            temp_state.transition[2] = 2 * mu_t
            temp_state.ss_prob = four_state_n.states[4].ss_prob
        three_state_n.states[j] = temp_state
    # Creates a four state system for adverse weather.
    four_state_a = System('4 State Transmission, Adverse')
    for i in range(1, 5):
        temp_state = State(i)
        if i == 1:
            temp_state.transition[2] = lambda_ta
            temp_state.transition[3] = lambda_ta
            temp_state.ss_prob = np.square(p_tua)
        elif i == 4:
            temp_state.transition[2] = mu_t
            temp_state.transition[3] = mu_t
            temp_state.ss_prob = np.square(p_tda)
        else:
            temp_state.transition[1] = mu_t
            temp_state.transition[4] = lambda_ta
            temp_state.ss_prob = p_tua * p_tda
        four_state_a.states[i] = temp_state
    # Creates a three state system for adverse weather.
    three_state_a = System('3 State Transmission, Adverse')
    for j in range(1, 4):
        temp_state = State(j)
        if j == 1:
            temp_state.transition[2] = 2 * lambda_ta
            temp_state.ss_prob = four_state_a.states[1].ss_prob
        elif j == 2:
            temp_state.transition[3] = lambda_ta
            temp_state.transition[1] = mu_t
            temp_state.ss_prob = four_state_a.states[2].ss_prob + four_state_a.states[3].ss_prob
        else:
            temp_state.transition[2] = 2 * mu_t
            temp_state.ss_prob = four_state_a.states[4].ss_prob
        three_state_a.states[j] = temp_state
    # Creates a 6 state system using the information from the two three state systems.
    six_state = System('6 State Transmission')
    for i in range(1, 7):
        temp_state = State(i)
        # Converts the 'normal' weather conditions to our first 3 states, and adds in the transition
        # to the adverse weather condition.
        if i in [1, 2, 3]:
            temp_state = three_state_n.states[i]
            temp_state.transition[i + 3] = n
        # Converts the 'adverse' weather conditions to our last 3 states, and adds in the transition
        # to the normal weather conditions.
        elif i in [4, 5, 6]:
            temp_state = three_state_a.states[i - 3]
            # Note that because the adverse weather transition rates were from states 1, 2, and 3, we need
            # to change them to represent the new state numbers, 4, 5, and 6. These are updated by creating
            # a new dictionary with the new state values, and replacing the old dictionary.
            temp_transition = dict()
            num_transition = len(temp_state.transition)
            for j in range(num_transition):
                curr_key = list(temp_state.transition.keys())[j]
                new_key = curr_key + 3
                rate_value = list(temp_state.transition.values())[j]
                temp_transition[new_key] = rate_value
            temp_state.transition = temp_transition
            temp_state.transition[i - 3] = s
        six_state.states[i] = temp_state
    # However, we need to update the probabilities of each state.
    # We do this through the BP = C method, where we get a transition rate matrix, and calculate the
    # steady state probabilities using that matrix.
    r_matrix = transition_rate(six_state)
    p_matrix = bpc(r_matrix)
    # Now we update the steady state probabiltiies for each state.
    for i in range(len(p_matrix)):
        six_state.states[i + 1].ss_prob = p_matrix[i]
    # Finally, we compress this system into a 3 state system, sorted by transmission capacity.
    three_state = System('3 State Transmission')
    for i in range(1, 4):
        temp_state = State(i)
        # Sums up the component states to create the new steady state probability.
        temp_state.ss_prob = six_state.states[i].ss_prob + six_state.states[i + 3].ss_prob
        # We know that in this instance, the component states are pairs, separated by 3 states.
        comp_1 = six_state.states[i]
        comp_2 = six_state.states[i + 3]
        # Hard coding calculation of new transition rates. Might revisit in the event I find a better way to do this.
        if i == 1:
            lambda_num = (comp_1.ss_prob * comp_1.transition[i + 1]) + \
                         (comp_2.ss_prob * comp_2.transition[i + 4])
            temp_state.transition[2] = lambda_num / temp_state.ss_prob
        elif i == 2:
            lambda_num = (comp_1.ss_prob * comp_1.transition[i + 1]) + \
                         (comp_2.ss_prob * comp_2.transition[i + 4])
            temp_state.transition[3] = lambda_num / temp_state.ss_prob
            mu_num = (comp_1.ss_prob * comp_1.transition[i - 1]) + \
                     (comp_2.ss_prob * comp_2.transition[i + 2])
            temp_state.transition[1] = mu_num / temp_state.ss_prob
        else:
            mu_num = (comp_1.ss_prob * comp_1.transition[i - 1]) + \
                     (comp_2.ss_prob * comp_2.transition[i + 2])
            temp_state.transition[2] = mu_num / temp_state.ss_prob
        three_state.states[i] = temp_state
    # Returns the final three state transmission system.
    return three_state


# Gets the steady state probabilities and transition rates for our final load states.
def load_info():
    lambda_4 = 6
    lambda_8 = 3
    # Creates a 5 state load system.
    load_state = System('5 State Load')
    for i in range(1, 6):
        temp_state = State(i)
        # State 5 is the only one that has a different transition rate, so we just set an if statement to catch it.
        if i == 5:
            temp_state.transition[1] = lambda_8
        # Otherwise, we set all the transition rates to be the same value.
        else:
            temp_state.transition[i + 1] = lambda_4
        load_state.states[i] = temp_state
    # In order to get the steady state probabilities, we calculate them using the BP = C method.
    r_matrix = transition_rate(load_state)
    p_matrix = bpc(r_matrix)
    # Update the steady state probabilities with the new values.
    for i in range(len(p_matrix)):
        load_state.states[i + 1].ss_prob = p_matrix[i]
    # Returns the final five state load system.
    return load_state


# Creates transition rate matrices.
def transition_rate(system):
    # Gets the number of states in the system, and initializes a zero matrix for the transition rate matrix.
    n = len(system.states)
    r_matrix = np.zeros((n, n))
    # Goes through each state in the system.
    for i in system.states:
        temp_state = system.states[i]
        # Gets the transition rates from that state, and adds it to to the r_matrix.
        for j in temp_state.transition:
            r_matrix[i - 1, j - 1] = temp_state.transition[j]
        # Constructs the diagonal terms by taking the negative of the sum of the off-diagonal terms in that row.
        r_matrix[i - 1, i - 1] = -sum(r_matrix[i - 1])
    # Returns the transition rate matrix.
    return r_matrix


# Calculates the steady state probabilities using the transition rate matrix.
def bpc(r_matrix):
    # Creates the b matrix, which is the transpose of the r matrix.
    b_mat = np.copy(np.transpose(r_matrix))
    # We set the first row to be all 1's, for the condition that the sum of the probabilities must be 1.
    # Technically, this can be done for any row, but for sake of simplicity, we'll focus only on the first row.
    b_mat[0, :] = 1
    # Create a C matrix, where all the elements are 0 except the first row, which is 1/
    c_mat = np.zeros(len(b_mat))
    c_mat[0] = 1
    # Solve the equation to get the steady state probabilities for each of our states.
    p_mat = np.linalg.solve(b_mat, c_mat)
    # Returns the steady state probability matrix.
    return p_mat
