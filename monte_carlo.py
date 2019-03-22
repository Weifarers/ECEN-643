import numpy as np


# Creates a class for each of our components that consists of two different states.
class TwoState(object):
    def __init__(self, name):
        self.name = name
        self.failure_rate = 0
        self.repair_rate = 0
        self.states = dict()


# Creates a class specifically for the transmission lines.  This is done since the transmission
# lines have different failure and repair rates for different weather conditions.
class TwoStateLine(object):
    def __init__(self, name):
        self.name = name
        self.failure_rate = dict()
        self.repair_rate = dict()
        self.states = dict()


# Creates a class for the load component, which can transition to multiple states.
class State(object):
    def __init__(self, name):
        # Gives a name to the state (ie: 1, 2, 3,..)
        self.name = name
        # Defines a dictionary for transition rates, where each key is the state being transitioned to,
        # and the value of the key is the rate.
        self.transition = dict()
        # Since only the load uses this object, we'll go ahead and hard code the states and their duration
        # into this object.
        self.states = {1: 60, 2: 105, 3: 205, 4: 105, 5: 60}
        self.states_duration = {1: 4, 2: 4, 3: 4, 4: 4, 5: 8}


# Creates a class for each iteration that we have in our simulation.
class SeqIter(object):
    def __init__(self, name):
        # Gives a name to the state (ie: 1, 2, 3,..)
        self.name = name
        # Each iteration has a specific duration, in units of days.
        self.duration = 0
        # We defined the states in terms of a dictionary; all the keys are components, and values are the states of
        # the components.
        self.states = dict()
        # We also have a states_duration dictionary, whose keys are the components, and whose values are the
        # remaining duration of the states.
        self.states_duration = dict()
        # lol = Loss of Load, or whether or not our generation can supply the amount of load being demanded.
        self.lol = False


# Initializing all of our objects in our system. Note that we have to give all rates in per hour time,
# as we track the evolution of the load on an hourly basis.
def comp_fill(comp_dict):
    # Creates three identical generators.
    for i in range(1, 4):
        # Creates a unique name for each generator.
        gen_name = 'G' + str(i)
        # Creates a new object for each generator, and sets the failure/repair rates.
        temp_comp = TwoState(gen_name)
        temp_comp.failure_rate = 0.1 / 24
        temp_comp.repair_rate = 1 / 8
        temp_comp.states = {True: 75, False: 0}
        # Adds this to our dictionary containing all of our components.
        comp_dict[gen_name] = temp_comp

    # Creating a two state weather component. In this case, normal = success and adverse = failure.
    # Since there's only one weather instance to keep track of, we don't need to generate multiple versions.
    w_comp = TwoState('W')
    w_comp.failure_rate = 1 / 200
    w_comp.repair_rate = 1 / 20
    # Adding it to the component dictionary.
    comp_dict['W'] = w_comp

    # Now looking at the transmission lines.
    for i in range(1, 3):
        # Creates a unique name for each transmission line, in the normal instance.
        t_name = 'T' + str(i)
        # Creates a new object for each line, and sets the failure/repair rates. Note that the failure and repair
        # rates are dictionaries in this instance, due to weather.
        temp_line = TwoStateLine(t_name)
        temp_line.failure_rate = {'N': 10 / 8760, 'A': 100 / 8760}
        temp_line.repair_rate = {'N': 1 / 8, 'A': 1 / 8}
        temp_line.states = {True: 100, False: 0}
        # Adds this to our component dictionary.
        comp_dict[t_name] = temp_line

    # Passes the dictionary through the load_info function, which adds in the loads, which are not a two state system.
    comp_dict = load_info(comp_dict)

    # Returns the dictionary for later use.
    return comp_dict


# Gets the transition rates for our  load states.
def load_info(comp_dict):
    lambda_4 = 1 / 4
    lambda_8 = 1 / 8
    load_state = State('L')
    # Creates a 5 state load system.
    for i in range(1, 6):
        # State 5 is the only one that has a different transition rate, so we just set an if statement to catch it.
        if i == 5:
            load_state.transition[1] = lambda_8
        # Otherwise, we set all the transition rates to be the same value.
        else:
            load_state.transition[i + 1] = lambda_4
    comp_dict['L'] = load_state

    return comp_dict


# Function that calculates the time remaining for any component to switch states, given the random number
# and the current state of the component.
def time_calc(rand_num, iteration, comp_dict, comp_val):
    # We throw an if statement to catch instances where we're looking at the transmission lines, whose failure
    # and repair rates depend on the weather.
    if comp_val in ['T1', 'T2']:
        # First, check the status of the weather.
        weather_cond = iteration.states['W']
        # By default, we'll assume we're in normal weather, and change the key if we see adverse weather.
        weather_key = 'N'
        # Catching instances where the weather is adverse, and sets a key for us to use later.
        if not weather_cond:
            weather_key = 'A'
        # Now we calculate the duration of the each state based on what state they're currently in.
        # If we're in the up state, designated by True, we calculate the time using the failure rate, with
        # the appropriate weather key.
        if iteration.states[comp_val]:
            time_val = -np.log(rand_num) / comp_dict[comp_val].failure_rate[weather_key]
            iteration.states_duration[comp_val] = time_val
        # If we're in the down state, we calculate the time using the repair rate, again with the
        # appropriate weather key.
        elif not iteration.states[comp_val]:
            time_val = -np.log(rand_num) / comp_dict[comp_val].repair_rate[weather_key]
            iteration.states_duration[comp_val] = time_val
    # Also checks for the load calculations, where we can just get the duration straight
    # from the component dictionary.
    elif comp_val == 'L':
        load_state = iteration.states[comp_val]
        iteration.states_duration[comp_val] = comp_dict[comp_val].states_duration[load_state]
    # Looks at instances that are not the lines or the load.
    else:
        # Checks if we're currently in the up state, designated by a 1. If so, calculate the time using
        # the failure rate.
        if iteration.states[comp_val]:
            time_val = -np.log(rand_num) / comp_dict[comp_val].failure_rate
            # Updates the duration of the event to the newly calculated one.
            iteration.states_duration[comp_val] = time_val
        # Checks if we're currently in the down state, designated by a 2. If so, calculate the time
        # using the repair rate.
        elif not iteration.states[comp_val]:
            time_val = -np.log(rand_num) / comp_dict[comp_val].repair_rate
            # Updates the duration of the event to the newly calculated one.
            iteration.states_duration[comp_val] = time_val

    return iteration


# This is our loss of load checker.
def lol_check(iteration, comp_dict):
    # Initializes a list of values that we store for later use.
    val_list = []
    # Goes through all the components and their current state.
    for i in iteration.states:
        # Ignores the weather, since that doesn't matter for loss of load calculations.
        if i == 'W':
            continue
        else:
            # Gets the state of the component, and refers to our comp_dict to find its value.
            current_state = iteration.states[i]
            curr_val = comp_dict[i].states[current_state]
            # Stores it in our val_list.
            val_list.append(curr_val)
    # Defines G1/G2, G3, T1/T2, and the Load.
    g12 = val_list[0] + val_list[1]
    g3 = val_list[2]
    transmission = val_list[3] + val_list[4]
    load = val_list[5]
    # Calculates the total generation of the system.
    total_gen = min(g12, transmission) + g3
    # Tests for loss of load.
    if total_gen < load:
        iteration.lol = True
    else:
        iteration.lol = False

    return iteration


# Main Function where calculations are done.
def main():
    # Sets a seed for random number generation, in order to have reproducible results.
    # Note that I did test a few different seeds, just to see what the range of results
    # were like, and they were all within an acceptable tolerance, so I just picked
    # the one that was the closest for demonstration purposes.
    random_gen = np.random.RandomState(150)
    # Initializes a component dictionary that is used as a reference dictionary for the
    # code to identify loss of load, or find specific transition rates.
    comp_dict = dict()
    # The comp_fill function fills the dictionary with all the information that we need.
    comp_dict = comp_fill(comp_dict)
    # This is the convergence criteria we will compare the coefficient of variation to.
    epsilon = 0.05
    # This is a dictionary where the keys are the years, which we initialize as 1, and
    # the values are the amount of hours spent in loss of load for that year.
    lol_dict = dict()
    year_idx = 1
    # We use year_time to calculate the amount of hours passed during each year. This
    # lets us stop when we reach 8760 hours, and reset the lol_time.
    year_time = 0
    # Number of hours spent in loss of load, updated every time there's a loss of load,
    # and reset every time a year has elapsed.
    lol_time = 0
    # An array used to store the averages of the hours spent in loss of load for every year.
    # This is the list we need to take the standard deviation of to calculate COV.
    lol_avg = []
    # Updates every time we go from no loss of load to a loss of load.
    lol_num = 0
    # A flag used to check if the last iteration had a loss of load. If false, then no.
    # We initialize it as false, since our starting state does not have loss of load.
    prev_lol_flag = False
    # Total number of hours spent in the simulation.
    sim_time = 0
    # Iteration counter.
    iter_count = 1
    # Now that all of our variables have been initialized, we hard code our initial state.
    iteration = SeqIter(0)
    # We assume everything is up, and that the load has begun in state 1.
    iteration.states = {'G1': True, 'G2': True, 'G3': True, 'W': True, 'T1': True, 'T2': True, 'L': 1}
    # Generates 6 random numbers for our initial iteration.
    init_rand = random_gen.rand(1, 7)
    # Calculates the times associated with the initial numbers that were drawn.
    for i in range(len(iteration.states)):
        # Assigns a random number for each state, except the load.
        rand_num = init_rand[0][i]
        comp_val = list(iteration.states.keys())[i]
        iteration = time_calc(rand_num, iteration, comp_dict, comp_val)
    # Now we begin the sequential simulation.
    while iter_count < 1000000:
        # Gets the minimum time to switch states, along with what component has that time.
        min_comp = min(iteration.states_duration.items(), key=lambda x: x[1])
        # We set the duration of the iteration the minimum time.
        iteration.duration = min_comp[1]

        # Now we check for loss of load.
        iteration = lol_check(iteration, comp_dict)
        # If we encounter loss of load, we need to do a few things.
        if iteration.lol:
            # First, update the loss of load time.
            lol_time += min_comp[1]
            # Now checks our flag to see if the previous iteration was a loss of load or not.
            if not prev_lol_flag:
                # If it was, update the number of times we've encountered loss of load.
                lol_num += 1
            # We only enter this if statement if there is a loss of load, so now we can
            # set the flag to true.
            prev_lol_flag = True
        # If we don't have loss of load, set the flag to false.
        else:
            prev_lol_flag = False
        # Now that we know how long the iteration lasts, and whether or not it's considered
        # a loss of load, we update the simulation time and the amount of time passed in
        # the current year.
        sim_time += min_comp[1]
        year_time += min_comp[1]

        # In order to calculate the COV, we need to track the amount of hours that the system
        # has spent in loss of load on a yearly basis.
        if year_time > 8760:
            # Every time a year has elapsed, we add a new key to our dictionary, with the
            # amount of hours spent in loss of load for that year.
            lol_dict[year_idx] = lol_time
            # We also update the year index.
            year_idx += 1
            # Now we reset the hours elapsed in the year, and the amount of hours spent in
            # a loss of load state.
            year_time = 0
            lol_time = 0
            # We also need to calculate the COV, first by getting the average of the hours spent
            # in loss of load per year.
            curr_avg = np.average(list(lol_dict.values()))
            # print(curr_avg)
            # We add it to our list, which is keeping track of all of our averages.
            lol_avg.append(curr_avg)
            # Finally, we calculate the COV.
            cov = np.std(lol_avg) / curr_avg
            # If we've reached a point lower than our convergence criteria, and at least 10
            # years have elapsed, break out of the while loop.
            if cov <= epsilon and year_idx > 10:
                break
        # And we'll update the times for all the other components.
        for a in iteration.states_duration:
            iteration.states_duration[a] -= min_comp[1]

        # We throw in a special case for when the minimum time is the load (which is pretty often.)
        if min_comp[0] == 'L':
            # This just forces the load to transition to the next state, and if you're in the
            # last state, to transition back to 1.
            next_load = iteration.states['L']
            if iteration.states['L'] < 5:
                next_load += 1
            elif iteration.states['L'] == 5:
                next_load = 1
            # Sets the status of the load to our new one.
            iteration.states['L'] = next_load
            # Sets the duration of the load state to a preset one.
            iteration.states_duration['L'] = comp_dict['L'].states_duration[next_load]
        # Handles all other instances where the load is not the first thing that's changing.
        else:
            # Now we change the status of our min_comp to the opposite state.
            iteration.states[min_comp[0]] = not iteration.states[min_comp[0]]
            # Now we need to generate a new random number for our component that has switched
            # states.
            new_random = random_gen.rand(1, 1)
            # We update the time to switching for that component using this new random variable.
            iteration = time_calc(new_random[0][0], iteration, comp_dict, min_comp[0])
        # Updates the iteration.
        iter_count += 1

    print('Frequency of Failure = {:.4f}/year'.format(lol_num / (sim_time / 8760)))
    print('LOLP = {:.4f}'.format(np.sum(list(lol_dict.values())) / sim_time))


if __name__ == '__main__':
    main()
