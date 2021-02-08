class Problem:
    def __init__(self,nVar=30,VarMin=-100,VarMax=100,minimize=1):
        self.nVar = nVar  # Number of Variables
        self.VarMin = VarMin  # Lower Bound of Decision Variables
        self.VarMax = VarMax  # Upper Bound of Decision Variables
        self.minimize = minimize  # Minimization 1 else Maximization 0


class Params:
    def __init__(self,problem,MaxIt=2000,MaxInteractions=2000,nPop=100,PopSize=100,w=2,wdamp=0.9,wmin=0.1,c1=2,c2=2,
                 ShowIterInfo=1,SelMethod=1,CrossMethod=1,CrossPercent=60,MutatPercent=20,EarlyExit=0):
        self.MaxIt = MaxIt  # PSO Max Interactions ( 0 if GA )
        self.MaxInteractions = MaxInteractions  # GA Max Interactions ( 0 if PSO )
        self.nPop = nPop  # PSO Population
        self.PopSize = PopSize  # GA Population
        self.w = w  # Inertia Coefficient
        self.wdamp = wdamp  # Damping Ratio of Inertia Coefficient
        self.wmin = wmin    # Minimum inertia value for linear decrease
        self.c1 = c1  # Personal Acceleration Coefficient
        self.c2 = c2  # Social Acceleration Coefficient
        self.ShowIterInfo = ShowIterInfo  # Flag for Showing Iteration Information
        self.SelMethod = SelMethod  # Breeding Method
        self.CrossMethod = CrossMethod  # Cross Breeding Method
        self.CrossPercent = CrossPercent  # Cross Breeding Percentage
        self.MutatPercent = MutatPercent  # Mutation Percentage
        self.ElitPercent = 100 - self.CrossPercent - self.MutatPercent  # Elite Population Percentage
        self.CrossNum = round(self.CrossPercent / 100 * self.PopSize)
        self.EarlyExit = EarlyExit  # Early exit criteria in effect ?

        if self.CrossNum % 2 != 0:
            self.CrossNum = self.CrossNum - 1

        self.MutatNum = round(self.MutatPercent / 100 * self.PopSize)
        self.ElitNum = self.PopSize - self.CrossNum - self.MutatNum

        if self.ElitNum % 2 != 0:
            self.ElitNum = self.ElitNum - 1

        self.MaxVelocity = 0.5 * (problem.VarMax - problem.VarMin)
        self.MinVelocity = -self.MaxVelocity
