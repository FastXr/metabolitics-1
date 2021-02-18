from cobra.flux_analysis import flux_variability_analysis
from cobra.util import fix_objective_as_constraint
from utils import load_network_model
from sympy.core.singleton import S
from .drugAnalysis import DrugReactionAnalysis
from cobra import Reaction
import extensions

class MetaboliticsAnalysis:
    '''
    Metabolitics analysis for metabolic dataset.
    '''

    def __init__(self, model='recon3D', drug='', target='', without_transports=True):
        '''
        :param model: cobra Model
        '''

        self.model = load_network_model(model=model)
        self.model.solver = 'cplex'
        # self.model.solver.configuration.timeout = 10 * 60
        self.without_transports = without_transports
        self.drug = drug
        self.target = target
        self.drug_analyzer = DrugReactionAnalysis()

    def drug_knock_out(self):
        for gene in self.drug_analyzer.drug_target(self.drug):
            try:
                self.model.genes.get_by_id(gene).knock_out()
            except:continue

    def combined_drug_knockout(self):
        # This is not needed
        for drug in self.drug:
            for gene in self.drug_analyzer.drug_target(drug):
                try:
                    self.model.genes.get_by_id(gene).knock_out()
                except:continue

    def set_objective(self, measured_metabolites):
        '''
        Updates objective function for given measured metabolites.

        :param dict measured_metabolites: dict in which keys are metabolite names 
            and values are float numbers represent fold changes in metabolites. 
        '''
        self.clean_objective()

        for k, v in measured_metabolites.items():

            m = self.model.metabolites.get_by_id(k)
            total_stoichiometry = m.total_stoichiometry(self.without_transports)
            for r in m.producers(self.without_transports):
                update_rate = v * r.metabolites[m] / total_stoichiometry
                r.objective_coefficient += update_rate


    def add_constraint(self, measured_metabolites):
        '''
        Add measurements as constraint to model.

        :param dict measured_metabolites: dict in which keys are metabolite names 
            and values are float numbers represent fold changes in metabolites.
        '''
        self.set_objective(measured_metabolites)
        fix_objective_as_constraint(self.model)


    def variability_analysis(self, measured_metabolites):
        if type(self.drug) == list:
            self.combined_drug_knockout()
        elif self.drug != '':
            self.drug_knock_out()
        elif self.target != '':
            self.model.genes.get_by_id(self.target).knock_out()

        self.set_objective(measured_metabolites)
        return flux_variability_analysis(self.model)

    def clean_objective(self):
        '''
        Cleans previous objective.
        '''
        self.model.objective = S.Zero

    def copy(self):
        return self.__class__(model=self.model.copy(), drug=self.drug, without_transports=self.without_transports)
