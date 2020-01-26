from typing import List

from guacamol_local.distribution_learning_benchmark import DistributionLearningBenchmark, ValidityBenchmark, \
    UniquenessBenchmark
from guacamol_local.goal_directed_benchmark import GoalDirectedBenchmark
from guacamol_local.scoring_function import ArithmeticMeanScoringFunction
from guacamol_local.standard_benchmarks import * 


def goal_directed_benchmark_suite(version_name: str) -> List[GoalDirectedBenchmark]:
    if version_name == 'v1':
        return goal_directed_suite_v1()
    if version_name == 'v2':
        return goal_directed_suite_v2()
    if version_name == 'v3':
        return goal_directed_suite_v3()
    if version_name == 'v4':
        return goal_directed_suite_v4()
    if version_name == 'v5':
        return goal_directed_suite_v5()
    if version_name == 'trivial':
        return goal_directed_suite_trivial()
    if version_name == 'test':
        return goal_directed_suite_test()
    if version_name == 'wsa':
        return goal_directed_suite_wsa()
    if version_name == 'wsc':
        return goal_directed_suite_wsc()
    if version_name == 'wosa':
        return goal_directed_suite_wosa()

    raise Exception(f'Goal-directed benchmark suite "{version_name}" does not exist.')


def distribution_learning_benchmark_suite(chembl_file_path: str,
                                          version_name: str,
                                          number_samples: int) -> List[DistributionLearningBenchmark]:
    """
    Returns a suite of benchmarks for a specified benchmark version

    Args:
        chembl_file_path: path to ChEMBL training set, necessary for some benchmarks
        version_name: benchmark version

    Returns:
        List of benchmaks
    """

    # For distribution-learning, v1 and v2 are identical
    if version_name == 'v1' or version_name == 'v2':
        return distribution_learning_suite_v1(chembl_file_path=chembl_file_path, number_samples=number_samples)

    raise Exception(f'Distribution-learning benchmark suite "{version_name}" does not exist.')


def goal_directed_suite_v1() -> List[GoalDirectedBenchmark]:
    max_logP = 6.35584
    return [
        isomers_c11h24(mean_function='arithmetic'),
        isomers_c7h8n2o2(mean_function='arithmetic'),
        isomers_c9h10n2o2pf2cl(mean_function='arithmetic', n_samples=100),

        hard_cobimetinib(max_logP=max_logP),
        hard_osimertinib(ArithmeticMeanScoringFunction),
        hard_fexofenadine(ArithmeticMeanScoringFunction),
        weird_physchem(),

        # start pop benchmark
        # e.g.
        start_pop_ranolazine(),

        # similarity Benchmarks

        # explicit rediscovery
        similarity(smiles='CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F', name='Celecoxib',
                   fp_type='ECFP4', threshold=1.0, rediscovery=True),
        similarity(smiles='Cc1c(C)c2OC(C)(COc3ccc(CC4SC(=O)NC4=O)cc3)CCc2c(C)c1O', name='Troglitazone',
                   fp_type='ECFP4', threshold=1.0, rediscovery=True),
        similarity(smiles='CN(C)S(=O)(=O)c1ccc2Sc3ccccc3C(=CCCN4CCN(C)CC4)c2c1', name='Thiothixene',
                   fp_type='ECFP4', threshold=1.0, rediscovery=True),

        # generate similar stuff
        similarity(smiles='Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl',
                   name='Aripiprazole', fp_type='FCFP4', threshold=0.75),
        similarity(smiles='CC(C)(C)NCC(O)c1ccc(O)c(CO)c1', name='Albuterol',
                   fp_type='FCFP4', threshold=0.75),
        similarity(smiles='COc1ccc2[C@H]3CC[C@@]4(C)[C@@H](CC[C@@]4(O)C#C)[C@@H]3CCc2c1', name='Mestranol',
                   fp_type='AP', threshold=0.75),

        logP_benchmark(target=-1.0),
        logP_benchmark(target=8.0),
        tpsa_benchmark(target=150.0),

        cns_mpo(max_logP=max_logP),
        qed_benchmark(),
        median_camphor_menthol(ArithmeticMeanScoringFunction)
    ]


def goal_directed_suite_v2() -> List[GoalDirectedBenchmark]:
    return [
        # explicit rediscovery
        similarity(smiles='CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F', name='Celecoxib', fp_type='ECFP4',
                   threshold=1.0, rediscovery=True),
        similarity(smiles='Cc1c(C)c2OC(C)(COc3ccc(CC4SC(=O)NC4=O)cc3)CCc2c(C)c1O', name='Troglitazone', fp_type='ECFP4',
                   threshold=1.0, rediscovery=True),
        similarity(smiles='CN(C)S(=O)(=O)c1ccc2Sc3ccccc3C(=CCCN4CCN(C)CC4)c2c1', name='Thiothixene', fp_type='ECFP4',
                   threshold=1.0, rediscovery=True),

        # generate similar stuff
        similarity(smiles='Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl', name='Aripiprazole', fp_type='ECFP4',
                   threshold=0.75),
        similarity(smiles='CC(C)(C)NCC(O)c1ccc(O)c(CO)c1', name='Albuterol', fp_type='FCFP4', threshold=0.75),
        similarity(smiles='COc1ccc2[C@H]3CC[C@@]4(C)[C@@H](CC[C@@]4(O)C#C)[C@@H]3CCc2c1', name='Mestranol',
                   fp_type='AP', threshold=0.75),

        # isomers
        isomers_c11h24(),
        isomers_c9h10n2o2pf2cl(),

        # median molecules
        median_camphor_menthol(),
        median_tadalafil_sildenafil(),

        # all other MPOs
        hard_osimertinib(),
        hard_fexofenadine(),
        ranolazine_mpo(),
        perindopril_rings(),
        amlodipine_rings(),
        sitagliptin_replacement(),
        zaleplon_with_other_formula(),
        valsartan_smarts(),
        decoration_hop(),
        scaffold_hop(),
    ]


def goal_directed_suite_v3() -> List[GoalDirectedBenchmark]:
    return [
        sa_logp_target(target=8),
        sa_hard_osimertinib(),
        sa_hard_fexofenadine(),
        sa_ranolazine_mpo(),
        sa_perindopril_rings(),
        sa_amlodipine_rings(),
        sa_sitagliptin_replacement(),
        sa_zaleplon_with_other_formula(),
        sa_valsartan_smarts(),
        sa_decoration_hop(),
        sa_scaffold_hop()
    ]

def goal_directed_suite_v4() -> List[GoalDirectedBenchmark]:
    return [
        hard_osimertinib(),
        hard_fexofenadine(),
        ranolazine_mpo(),
        perindopril_rings(),
        amlodipine_rings(),
        sitagliptin_replacement(),
        zaleplon_with_other_formula(),
        valsartan_smarts(),
        decoration_hop(),
        scaffold_hop(),
    ]

def goal_directed_suite_v5() -> List[GoalDirectedBenchmark]:
    return [
        sc_hard_osimertinib(),
        sc_hard_fexofenadine(),
        sc_ranolazine_mpo(),
        sc_perindopril_rings(),
        sc_amlodipine_rings(),
        sc_sitagliptin_replacement(),
        sc_zaleplon_with_other_formula(),
        sc_valsartan_smarts(),
        sc_decoration_hop(),
        sc_scaffold_hop(),
    ]

def goal_directed_suite_test() -> List[GoalDirectedBenchmark]:
    return [
        qed_benchmark()
    ]


def goal_directed_suite_trivial() -> List[GoalDirectedBenchmark]:
    """
    Trivial goal-directed benchmarks from the paper.
    """
    return [
        logP_benchmark(target=-1.0),
        logP_benchmark(target=8.0),
        tpsa_benchmark(target=150.0),
        cns_mpo(),
        qed_benchmark(),
        isomers_c7h8n2o2(),
        pioglitazone_mpo(),
    ]

def goal_directed_suite_wosa() -> List[GoalDirectedBenchmark]:
    """
    Trivial goal-directed benchmarks from the paper.
    """
    return [
        cns_mpo(),
        qed_benchmark(),
        isomers_c7h8n2o2(),
        pioglitazone_mpo(),
    ]

def goal_directed_suite_wsa() -> List[GoalDirectedBenchmark]:
    """
    Trivial goal-directed benchmarks from the paper.
    """
    return [
        sa_cns_mpo(),
        sa_qed_benchmark(),
        sa_isomers_c7h8n2o2(),
        sa_pioglitazone_mpo(),
    ]

def goal_directed_suite_wsc() -> List[GoalDirectedBenchmark]:
    """
    Trivial goal-directed benchmarks from the paper.
    """
    return [
        sc_cns_mpo(),
        sc_qed_benchmark(),
        sc_isomers_c7h8n2o2(),
        sc_pioglitazone_mpo(),
    ]

def distribution_learning_suite_v1(chembl_file_path: str, number_samples: int = 10000) -> \
        List[DistributionLearningBenchmark]:
    """
    Suite of distribution learning benchmarks, v1.

    Args:
        chembl_file_path: path to the file with the reference ChEMBL molecules

    Returns:
        List of benchmarks, version 1
    """
    return [
        ValidityBenchmark(number_samples=number_samples),
        UniquenessBenchmark(number_samples=number_samples),
        novelty_benchmark(training_set_file=chembl_file_path, number_samples=number_samples),
        kldiv_benchmark(training_set_file=chembl_file_path, number_samples=number_samples),
        frechet_benchmark(training_set_file=chembl_file_path, number_samples=number_samples)
    ]
