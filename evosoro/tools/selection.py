import random
import math
import difflib


def pareto_selection(population):
    """Return a list of selected individuals from the population.

    All individuals in the population are ranked by their level, i.e. the number of solutions they are dominated by.
    Individuals are added to a list based on their ranking, best to worst, until the list size reaches the target
    population size (population.pop_size).

    Parameters
    ----------
    population : Population
        This provides the individuals for selection.

    Returns
    -------
    new_population : list
        A list of selected individuals.

    """
    new_population = []

    # SAM: moved this into calc_dominance()
    # population.sort(key="id", reverse=False) # <- if tied on all objectives, give preference to newer individual

    # (re)compute dominance for each individual
    population.calc_dominance()

    # sort the population multiple times by objective importance
    population.sort_by_objectives()

    # divide individuals into "pareto levels":
    # pareto level 0: individuals that are not dominated,
    # pareto level 1: individuals dominated one other individual, etc.
    done = False
    pareto_level = 0
    while not done:
        this_level = []
        size_left = population.pop_size - len(new_population)
        for ind in population:
            if len(ind.dominated_by) == pareto_level:
                this_level += [ind]

        # add best individuals to the new population.
        # add the best pareto levels first until it is not possible to fit them in the new_population
        if len(this_level) > 0:
            if size_left >= len(this_level):  # if whole pareto level can fit, add it
                new_population += this_level

            else:  # otherwise, select by sorted ranking within the level
                new_population += [this_level[0]]
                while len(new_population) < population.pop_size:
                    random_num = random.random()
                    log_level_length = math.log(len(this_level))
                    for i in range(1, len(this_level)):
                        if math.log(i) / log_level_length <= random_num < math.log(i + 1) / log_level_length and \
                                        this_level[i] not in new_population:
                            new_population += [this_level[i]]
                            continue

        pareto_level += 1
        if len(new_population) == population.pop_size:
            done = True

    for ind in population:
        if ind in new_population:
            ind.selected = 1
        else:
            ind.selected = 0

    return new_population


def pareto_tournament_selection(population):
    """Reduce the population pairwise.

    Two individuals from the population are randomly sampled and the inferior individual is removed from the population.
    This process repeats until the population size is reduced to either the target population size (population.pop_size)
    or the number of non-dominated individuals / Pareto front (population.non_dominated_size).

    Parameters
    ----------
    population : Population
        This provides the individuals for selection.

    Returns
    -------
    new_population : list
        A list of selected individuals.

    """
    # population.add_random_individual()  # adding in random ind in algorithms.py
    population.calc_dominance()
    random.shuffle(population.individuals)
    print "The nondominated size is", population.non_dominated_size

    while len(population) > population.pop_size and len(population) > population.non_dominated_size:

        inds = random.sample(range(len(population)), 2)
        ind0 = population[inds[0]]
        ind1 = population[inds[1]]

        if population.dominated_in_multiple_objectives(ind0, ind1):
            print "(fit) {0} dominated by {1}".format(ind0.fitness, ind1.fitness)
            print "(age) {0} dominated by {1}".format(ind0.age, ind1.age)
            population.pop(inds[0])
        elif population.dominated_in_multiple_objectives(ind1, ind0):
            print "(fit) {1} dominated by {0}".format(ind0.fitness, ind1.fitness)
            print "(age) {1} dominated by {0}".format(ind0.age, ind1.age)
            population.pop(inds[1])
        # else:
        #     population.pop(random.choice(inds))

    population.sort_by_objectives()

    return population.individuals

def pareto_selection_diversify_ancestry(population):
    """Return a list of selected individuals from the population.
        DANIEL: added comparison of lineage for each added individual. Before an inidividual is added, the lineage is
        compared to the already added population. If the similarities are higher than x%, the individual is not added
    All individuals in the population are ranked by their level, i.e. the number of solutions they are dominated by.
    Individuals are added to a list based on their ranking, best to worst, until the list size reaches the target
    population size (population.pop_size).

    Parameters
    ----------
    population : Population
        This provides the individuals for selection.

    Returns
    -------
    new_population : list
        A list of selected individuals.

    """
    new_population = []
    difference = []
    ratio = []
    # SAM: moved this into calc_dominance()
    # population.sort(key="id", reverse=False) # <- if tied on all objectives, give preference to newer individual

    # (re)compute dominance for each individual
    population.calc_dominance()

    # sort the population multiple times by objective importance
    population.sort_by_objectives()

    # divide individuals into "pareto levels":
    # pareto level 0: individuals that are not dominated,
    # pareto level 1: individuals dominated one other individual, etc.
    done = False
    pareto_level = 0
    rejected = 0
    while not done:
        this_level = []
        size_left = population.pop_size - len(new_population)
        for ind in population:
            if len(ind.dominated_by) == pareto_level:
                this_level += [ind]

        # add best individuals to the new population.
        # add the best pareto levels first until it is not possible to fit them in the new_population
        if len(this_level) > 0:
            if size_left >= len(this_level):  # if whole pareto level can fit, add it
                new_population += [this_level[0]]
                for i in range(0, len(this_level)):
                    if this_level[i] not in new_population:
                        for j in range(0, len(new_population)):
                            if len(population.lineage_dict) != 0:
                                difference.append(difflib.SequenceMatcher(None, population.lineage_dict.values()[i],
                                                                        population.lineage_dict.values()[j]))
                                ratiobuffer = difference[j].ratio()
                                if ratiobuffer==1:
                                    ratio.append(0)
                                else:
                                    ratio.append(ratiobuffer)
                            else:
                                ratio = [0,0]
                        if (sum(f>0.5 for f in ratio))/float(len(ratio))<0.25 or (sum(f>0.5 for f in ratio)) < 4: #if 0.5 of ancestors are the same for more then 0.1 of the individuals, don't select this individual
                            new_population += [this_level[i]]
                        else: rejected +=0.000000000001
                        continue
            else:  # otherwise, select by sorted ranking within the level
                new_population += [this_level[0]]
                while len(new_population) < population.pop_size and rejected<(len(this_level)-len(new_population)): #because individuals can be deselected, a pareto level can still be smaller then the size left after evaluation
                    #random_num = random.random()
                    #log_level_length = math.log(len(this_level))
                    for i in range(1, len(this_level)):
                        if this_level[i] not in new_population:
                            for j in range(0, len(new_population)):
                                if len(population.lineage_dict)!=0:
                                    difference.append(difflib.SequenceMatcher(None, population.lineage_dict.values()[i],population.lineage_dict.values()[j]))
                                    ratiobuffer = difference[j].ratio()
                                    if ratiobuffer == 1:
                                        ratio.append(0)
                                    else:
                                        ratio.append(ratiobuffer)
                                else: ratio=[0,0]
                            if (((sum(f>0.5 for f in ratio))/float(len(ratio)))<0.25 or (sum(f>0.5 for f in ratio)) < 4) and len(new_population) < population.pop_size:
                                new_population += [this_level[i]]
                            else: rejected +=1.0000000000
                            continue

        pareto_level += 1
        rejected =0
        if len(new_population) == population.pop_size:
            done = True

    for ind in population:
        if ind in new_population:
            ind.selected = 1
        else:
            ind.selected = 0

    return new_population

