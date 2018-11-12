import random
import math
import difflib
import numpy as np
import subprocess as sub

global fitness_list
global fitness_list_relative
global reset_counter
reset_counter =0
fitness_list= []
fitness_list_relative = []

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


def annealing_selection(population):
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
    max_gens = 1000 #hardcoded right now
    new_population = []
    # SAM: moved this into calc_dominance()
    # population.sort(key="id", reverse=False) # <- if tied on all objectives, give preference to newer individual

    # (re)compute dominance for each individual
    population.calc_dominance()

    # sort the population multiple times by objective importance
    population.sort_by_objectives()

    done = False
    annealing_temperature = (population.gen/max_gens +1.0)*6.0 -5.0
    #from first individual temp =1 to last individual temp=2
    population_buffer = []
    for ind in population:
        population_buffer += [ind] #to delete individuals that have been selected
    while not done:

        #using a beta distribution, for temp =1 it is a uniform distribution to
        # temp =2 a
        population_buffer_size = len(population_buffer)
        random_num = int(np.random.beta(1.0,annealing_temperature)*population_buffer_size)
        if population_buffer[random_num] not in new_population:
            new_population += [population_buffer[random_num]]
            del population_buffer[random_num]
        if len(new_population) == population.pop_size:
            done = True

    for ind in population:
        if ind in new_population:
            ind.selected = 1
        else:
            ind.selected = 0
    print "Beta is now ", annealing_temperature
    return new_population


def pareto_selection_reset(population):
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
    difference = 0.0
    ratio = 0.0

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
    population_buffer = []

    for ind in population:
        population_buffer += [ind] #to delete individuals that have been selected

    global reset_counter
    global fitness_list_relative
    global fitness_list
    deleted_individuals = []
    if reset_counter >=3: #301 counter after how many gens the first evaluation of local optima has to be done, when first reset is possible
        print "Sum of relative change in fitness over the last 300 generations is {0}".format(sum(fitness_list_relative[(len(fitness_list_relative) - 1):]))
        if sum(fitness_list_relative[(len(fitness_list_relative) - 1):])<0.0001: # sum of relative change in fitness over the last x generations
            j = population_buffer[0].parent_id #the parent of the best individual, to check the lineages as they are not updated yet
            if not j==-1: #tough chance, but the best individual is newly generated
                if population_buffer[0].id in population.lineage_dict:
                    j = population_buffer[0].id
                i=1
                while i < len(population_buffer):#for i in range(1, len(population_buffer)):
                    k = population_buffer[i].parent_id #get the parent id of the individual to be checked with the best individual
                    l = population_buffer[i].id
                    ratio = 0
                    if not k==-1: #check if the individual has more then one ancestor
                        if l in population.lineage_dict: #if the individual has a registred lineage use this one, otherwise the parents lineage
                            difference = difflib.SequenceMatcher(None, population.lineage_dict[l], population.lineage_dict[j])
                            ratio = difference.ratio()
                        elif not len(population.lineage_dict[k])==0: #check if the individuals parent has ancestors
                            difference = difflib.SequenceMatcher(None, population.lineage_dict[k], population.lineage_dict[j])
                            ratio = difference.ratio()  #else: ratio =0 # if lineage of parent =0
                    if ratio >= 0.8: # if there's a higher similarity then this number, the individual is close family of the best individual
                        deleted_individuals.append(l)
                        del population_buffer[i] # and is deleted  the population buffer shifted one place because of the deletion, so the same number has to be checked again., i stays the same, otherwise
                    else:
                        i += 1
                print "Deleted best individual {0} and close family {1}".format(population_buffer[0].id, deleted_individuals)
                deleted_individuals.append(population_buffer[0].id)
                run_directory = 'transmission_reset_data'
                run_name = 'transmission_reset'
                sub.call("mkdir " + run_directory + "/bestSoFar/LocalOptima/Gen_%04i" % population.gen, shell=True)
                for individual in population:
                    if individual.id in deleted_individuals:
                        sub.call(
                            "mv " + run_directory + "/ancestors/" + run_name + "--id_%05i.vxa" % individual.id +
                            " " + run_directory + "/bestSoFar/LocalOptima/Gen_%04i/" % population.gen + "/" +
                            run_name + "Gen_%04i--Fit_%.08f--id_%05i--dom_%d.vxa" %
                            (population.gen, individual.fitness, individual.id, len(individual.dominated_by)), shell=True)

                del population_buffer[0] #after every one is checked, the best individual is deleted
                individual = population_buffer[0] #add the new best result in the folder as well
                sub.call(
                    "mv " + run_directory + "/ancestors/" + run_name + "--id_%05i.vxa" % individual.id +
                    " " + run_directory + "/bestSoFar/LocalOptima/Gen_%04i/" % population.gen + "/" +
                    run_name + "Gen_%04i--Fit_%.08f--id_%05i--dom_%d_new_best_id.vxa" %
                    (population.gen, individual.fitness, individual.id, len(individual.dominated_by)), shell=True)
                reset_counter=0 #and the counter is reset to let it run to find a new local optimum
                fitness_list = []
                fitness_list_relative = []
                population.best_fit_so_far = 0
    else:
        reset_counter+=1


    while not done:
        population_buffer_size = len(population_buffer)
        this_level = []
        size_left = population.pop_size - len(new_population)
        for ind in population_buffer:
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
        elif pareto_level > population.pop_size:
            done = True
            print "Deleted to many individuals in reset to have full population."

    for ind in population:
        if ind in new_population:
            ind.selected = 1
        else:
            ind.selected = 0
    fitness_list += [0]
    fitness_list_relative += [0]
    fitness_list[population.gen] = population.best_fit_so_far #list of best individuals per generation
    fitness_list_relative[population.gen] = (population.best_fit_so_far-fitness_list[population.gen-1])/fitness_list[population.gen-1]
    return new_population