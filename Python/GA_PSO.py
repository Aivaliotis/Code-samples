from input_data import Problem, Params
import numpy as np
import time
# import cProfile


def cost1(x):
    cost = sum(x**2)
    return cost


def GA(problem, params, CostFunct):
    start_time = time.time()

    # Initializon of functions
    class InitialGlobalBest:
        position = []
        value = np.inf

    class Individual:
        def __init__(self):
            self.position = np.random.rand(problem.nVar) * (problem.VarMax - problem.VarMin) + problem.VarMin  # individual's position
            self.value = []  # value of individual

        # evaluate current fitness
        def evaluate(self):
            self.value = CostFunct(self.position)

    def mutation(par1):
        new_gen = par1
        new_gen.position = -par1.position
        return new_gen

    def cross_breed(par1, par2):
        children = []
        for i in range(2):
            children.append(Individual())
            if params.CrossMethod == 1:
                a = np.random.rand(problem.nVar)
            elif params.CrossMethod == 2:
                a = np.random.rand(1)

            children[i].position = (a * par1.position) + (1 - a) * par2.position

        return children

    # establish the Population
    individual = []
    GlobalBest = InitialGlobalBest()
    value = []

    for i in range(0, params.PopSize):
        individual.append(Individual())
        individual[i].evaluate()

        if individual[i].value < GlobalBest.value:
            GlobalBest.position = individual[i].position
            GlobalBest.value = individual[i].value

        value.append(individual[i].value)

    interaction = 0
    # main Loop
    while interaction < params.MaxInteractions:
        new_gen = individual
        ParentIndex = np.argsort(value)

        elite_pop = ParentIndex[range(0, int(params.ElitNum / 2))]
        for i in range(0, int(params.ElitNum / 2)):
            new_gen[2*i:2*i+1] = cross_breed(individual[elite_pop[i]], individual[elite_pop[-i]])

        cross_pop = ParentIndex[range(0, int(params.CrossNum / 2))]
        for i in range(0, int(params.CrossNum / 2)):
            new_gen[params.ElitNum+2*i:params.ElitNum+2*i + 1] = cross_breed(individual[cross_pop[i]], individual[cross_pop[-i]])

        for i in range(params.MutatNum):
            new_gen[params.ElitNum+params.CrossNum+i] = mutation(individual[np.random.randint(0,params.PopSize)])

        for i in range(0, params.PopSize):
            individual[i] = new_gen[i]
            individual[i].evaluate()

            if individual[i].value < GlobalBest.value:
                GlobalBest.position = individual[i].position
                GlobalBest.value = individual[i].value

            value[i] = individual[i].value

        interaction += 1
        if params.ShowIterInfo:
            print('GA: Iteration = ', interaction, ' GlobalBest = ', GlobalBest.value,' time =', time.time() - start_time, 's')

    return GlobalBest, time.time() - start_time

def PSO(problem, params, CostFunct):

    start_time = time.time()
    # Initialize Global Best
    class InitialGlobalBest:
        position = []
        value = np.inf

    class Particle:
        def __init__(self):
            self.position = np.random.rand(problem.nVar) * (problem.VarMax - problem.VarMin) + problem.VarMin  # particle position
            self.velocity = np.zeros(problem.nVar)  # particle velocity
            self.best_position = self.position  # best position individual
            self.best_value = np.inf  # best value of individual
            self.value = []  # value individual

        # evaluate current fitness
        def evaluate_and_update_Best(self, GlobalBest):
            self.value = CostFunct(self.position)

            # check to see if the current position is an individual best
            if self.value < self.best_value:
                self.best_position = self.position.copy()
                self.best_value = self.value
                if self.value < GlobalBest.value:
                    GlobalBest.position = self.position
                    GlobalBest.value = self.value

        # update new particle velocity
        def update_velocity_and_position(self, GlobalBest, it):
            rand1 = rand2 = []
            if params.CrossMethod == 1:
                rand1 = np.random.rand(problem.nVar)
                rand2 = np.random.rand(problem.nVar)
            elif params.CrossMethod == 2:
                rand1 = np.random.rand(1)
                rand2 = np.random.rand(1)

            # Update Velocity
            personal_velocity = params.c1 * rand1 * (self.best_position - self.position)
            social_velocity = params.c2 * rand2 * (GlobalBest.position - self.position)
            if params.wmin <= 0:
                self.velocity = (params.w - params.wmin) * (it/params.MaxIt) * self.velocity + personal_velocity + social_velocity
            else:
                self.velocity = params.w * (params.wdamp ** it) * self.velocity + personal_velocity + social_velocity

            # Apply Velocity boundary conditions
            self.velocity = np.minimum(self.velocity, params.MaxVelocity)
            self.velocity = np.maximum(self.velocity, params.MinVelocity)

            # Update Position
            self.position = self.position + self.velocity

            # Apply Position boundary conditions
            self.position = np.minimum(self.position, problem.VarMax)
            self.position = np.maximum(self.position, problem.VarMin)

    # establish the Population
    particle = []
    GlobalBest = InitialGlobalBest()
    for i in range(0, params.PopSize):
        particle.append(Particle())
        particle[i].evaluate_and_update_Best(GlobalBest)

    it = 0
    # main Loop
    while it < params.MaxIt:
        for i in range(0, params.nPop):
            particle[i].update_velocity_and_position(GlobalBest, it)
            particle[i].evaluate_and_update_Best(GlobalBest)

        it += 1

        if params.ShowIterInfo:
            print('PSO: Iteration = ', it, ' GlobalBest = ', GlobalBest.value, ' time =', time.time() - start_time, 's')

    return GlobalBest, time.time() - start_time


if __name__ =='__main__':
    problem = Problem()
    params = Params(problem)
    # pr = cProfile.Profile()
    # pr.enable()
    PSO_sol, PSO_time = PSO(problem, params, cost1)
    GA_sol, GA_time = GA(problem, params, cost1)

    print('PSO ', 'Value = ', PSO_sol.value, ' time =', PSO_time, 's')
    print('GA ', 'Value = ', GA_sol.value, ' time =', GA_time, 's')
    # pr.disable()
    # pr.print_stats()
