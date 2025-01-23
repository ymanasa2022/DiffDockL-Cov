######
# By Manasa Yadavalli 
# Jan 2025
######
import glob
import os 
from posebusters import PoseBusters
from pathlib import Path
import pandas as pd
import numpy as np 
from rdkit import Chem
from posebusters.modules.rmsd import check_rmsd

# use posebusters conda env
def posebusters_df(true_poses, pred_poses, pose_bust_cols):
    '''
    performs posebusters tests on predicted poses of all proteins 
    inputs
    -------
    true_poses (str): path to ground truth poses 
    pred_poses(str): path to predicted poses
    pose_bust_cols(list): column names of pose busters output to make final df

    outputs
    -------
    all_pose_evals: returns a DataFrame with the posebusters output for each predicted pose for each protein
    '''
    # import ipdb; ipdb.set_trace()
    all_pose_evals = pd.DataFrame(columns=pose_bust_cols)
    for prot in glob.glob(f"{true_poses}/**/*.pdb", recursive=True):
        prot_path_split = prot.split('/') # protein file path dir list
        complex_name = prot_path_split[-1].split('_')[0] # pdb code 
        print(complex_name)

        true_lig = f'{true_poses}/{complex_name}/{complex_name}_ligand.sdf'
        poses = glob.glob(f"{pred_poses}/{complex_name}/*.sdf")

        rank1_file = [pose for pose in poses if 'rank1.sdf' in pose]
        if len(poses) < 10: # to filter out failed sampling attempts 
            print(f'not all poses were predicted for {complex_name}, moving to next')
        else:
            poses.remove(rank1_file[0]) # removing duplicate rank1 pose
            dfs = []
            for pose in poses: 
                
                buster = PoseBusters(config="redock")
                buster_current = buster.bust(pose, true_lig, prot, full_report=True)

                df_buster_reset = buster_current.reset_index() # changing file and molecule to columns instead of indexes
                file = df_buster_reset.iloc[0, df_buster_reset.columns.get_loc('file')]
                rank_num = file.split('/')[-1].split('_')[0]
                # adding complex and pose rank to df 
                df_buster_reset.at[0, 'molecule'] = f'{complex_name}_{rank_num}'
                dfs.append(df_buster_reset)

            all_pose_evals = pd.concat([all_pose_evals] + dfs, ignore_index=True)
    return all_pose_evals

def rmsd_only(true_poses, pred_poses, rmsd_threshold=2.0):
    '''
    performs rmsd check on predicted top-1 poses generated by DiffDock-L's score model
    inputs
    -------
    true_poses (str): path to ground truth poses 
    pred_poses (str): path to predicted poses
    rmsd_threshold (float): RMSD cutoff. Default: 2.0

    outputs
    -------
    top1_rmsd_perc(float): % of proteins with top-1 pose < rmsd_threshold
    median_top1 (float): median RMSD of all protein top-1 poses 
    std_top1 (float): standard deviation of the 
    no_poses (lst): list of proteins that did not have sufficient generated poses 
    metrics_df (DataFrame): assessment metrics for each protein 
    '''
    # initializing values/df
    columns = ['Protein', 'Best RMSD', 'Mean RMSD', 'Std Dev of RMSD', 'RMSD < 2 %']
    metrics_df = pd.DataFrame(columns=columns)

    rmsd_top1 = np.array([]) 
    prots = np.array([])
    no_poses = []
    for prot in glob.glob(f"{true_poses}/**/*.pdb", recursive=True):
        prot_path_split = prot.split('/') 
        complex_name = prot_path_split[-1].split('_')[0] # pdb code 
        # print(complex_name)

        true_lig = f'{true_poses}/{complex_name}/{complex_name}_ligand.sdf'
        poses = glob.glob(f"{pred_poses}/{complex_name}/*.sdf")

        rank1_file = [pose for pose in poses if 'rank1.sdf' in pose]
        if len(poses) < 10: # to filter out failed sampling attempts 
            print(f'not all poses were predicted for {complex_name}, moving to next')
            no_poses.append(prot)
            continue # move to next prot 

        poses.remove(rank1_file[0]) # removing duplicate rank1 pose

        count_poses = 0
        pose_rmsd = np.array([])
        for pose in poses: 
            true_ligand_supplier = Chem.SDMolSupplier(true_lig)
            predicted_pose_supplier = Chem.SDMolSupplier(pose)

            mol_true = true_ligand_supplier[0] # mol obj
            mol_pred = predicted_pose_supplier[0]
            rmsd_dict = check_rmsd(mol_pred, mol_true, rmsd_threshold)
            rmsd_result = rmsd_dict['results']['rmsd']

            pose_rmsd = np.append(pose_rmsd, rmsd_result)
            
            count_poses += (rmsd_result < rmsd_threshold)

        # append best pose's rmsd 
        rmsd_top1 = np.append(rmsd_top1, np.min(pose_rmsd))
        prots = np.append(prots, prot)

        # update metrics per prot 
        metrics_per_prot = pd.concat([metrics_df, pd.DataFrame({'Protein': prot,
        'Best RMSD': rmsd_top1,
        'Mean RMSD': np.mean(pose_rmsd),
        'Std Dev of RMSD': np.std(pose_rmsd),
        'RMSD < 2 %': count_poses/len(poses) * 100 , })], ignore_index=True)

    # median of top-1 RMSDs
    median_top1 = round(np.mean(rmsd_top1), 2)
    std_top1 = round(np.std(rmsd_top1), 2)

    # Top-1 RMSD % < 2 
    count = np.sum(rmsd_top1 < 2)
    top1_rmsd_perc = round(count/len(rmsd_top1) * 100, 2)

    # fine best prot, worst prot 
    prot_rmsd_best = dict(zip(prots, rmsd_top1))
    prot_min_rmsd = min(prot_rmsd_best, key=prot_rmsd_best.get)
    prot_max_rmsd= max(prot_rmsd_best, key=prot_rmsd_best.get)
    
    print('prot with min rmsd:', prot_min_rmsd, '\nrmsd:', round(prot_rmsd_best[prot_min_rmsd], 3))
    print('prot with max rmsd:', prot_max_rmsd, '\nrmsd:', round(prot_rmsd_best[prot_max_rmsd], 3))

    return top1_rmsd_perc, median_top1, std_top1, no_poses, metrics_per_prot

def score_model_posebust(all_pose_buster_df, num_gen=10, rmsd_check=2):
    '''
    assesses the performance of poses generated by DiffDock-L's score model
    inputs
    -------
        all_pose_buster_df (DataFrame): output from make_eval_df made using PoseBusters 
        rmsd_check (float): Å RMSD threshold check. Default: 2
        num_gen (int): number of poses generated for a ligand. Default: 10
    outputs
    -------
        best_rmsd_perc (int): % of lowest rmsd poses that have rmsd less than threshold 
        mean_rmsd_perc (int): % of mean rmsd's of all generated poses that have rmsd less than threshold 
        metrics_df (DataFrame): assessment metrics for each protein 
    '''
    all_pose_buster_df['protein'] = all_pose_buster_df['molecule'].str.extract(r'^(.*?)_') 
    all_pose_buster_df['rank'] = all_pose_buster_df['molecule'].str.extract(r'rank(\d+)').astype(int)
    # sort within proteins based on rmsd of poses generated NOT by rank predicted
    sorted_df = all_pose_buster_df.groupby('protein', group_keys=False).apply(lambda group: group.sort_values('rmsd'))
    group_df = sorted_df.groupby('protein')

    # initializing values/df
    columns = ['Protein', 'Best RMSD', 'Mean RMSD', 'Std Dev of RMSD', 'RMSD < 2 %']
    metrics_df = pd.DataFrame(columns=columns)

    rmsd_best_np = np.array([])
    for prot, row_index in group_df.groups.items(): 
        # best rmsd: top-1
        rmsd_best = sorted_df.loc[row_index[0], 'rmsd']
        rmsd_best_np = np.append(rmsd_best_np, rmsd_best)

        # mean rmsd assessment for each prot
        group_rmsd = group_df.get_group(prot)
        mean_rmsd = group_rmsd['rmsd'].mean()
        # spread of RMSDs for each prot
        std_rmsd = group_rmsd['rmsd'].std()

        # % of poses that has < 2 RMSD per prot
        count_poses = 0
        for index in row_index: 
            rmsd_value = sorted_df.loc[index, 'rmsd']
            count_poses += (rmsd_value < rmsd_check)
        rmsd_per_prot = count_poses/num_gen * 100 

        # update metrics
        metrics_per_prot = pd.concat([metrics_df, pd.DataFrame({'Protein': [prot],
        'Best RMSD': [rmsd_best],
        'Mean RMSD': [mean_rmsd],
        'Std Dev of RMSD': [std_rmsd],
        'RMSD < 2 %': [rmsd_per_prot], })], ignore_index=True)

    # median of top-1 RMSDs
    median_top1 = round(np.mean(rmsd_best_np), 2)
    std_top1 = round(np.std(rmsd_best_np), 2)

    # Top-1 RMSD % < 2 
    count = np.sum(rmsd_best_np < 2)
    print(rmsd_best_np)
    top1_rmsd_perc = round(count/len(rmsd_best_np) * 100, 2)

    return top1_rmsd_perc, median_top1, std_top1, metrics_per_prot

def confidence_model_eval(all_pose_buster_df, num_gen=10, rmsd_check=2):
    '''
    assesses the performance of poses ranked by DiffDock-L's confidence model
    input: 
        all_pose_buster_df (DataFrame): output from make_eval_df made using PoseBusters 
        rmsd_check (float): Å RMSD threshold check. Default: 2
        num_gen (int): number of poses generated for a ligand
    output: 
        top_pose: % of times when the top ranked pose had the lowest rmsd amongst all sampled poses 
    '''
    # best pose 
    

    return rank1_corrects

def main():
    pose_bust_cols = ['file',
    'molecule',
    'mol_pred_loaded',
    'mol_true_loaded',
    'mol_cond_loaded',
    'sanitization',
    'inchi_convertible',
    'all_atoms_connected',
    'molecular_formula',
    'molecular_bonds',
    'double_bond_stereochemistry',
    'tetrahedral_chirality',
    'bond_lengths',
    'bond_angles',
    'internal_steric_clash',
    'aromatic_ring_flatness',
    'double_bond_flatness',
    'internal_energy',
    'protein-ligand_maximum_distance',
    'minimum_distance_to_protein',
    'minimum_distance_to_organic_cofactors',
    'minimum_distance_to_inorganic_cofactors',
    'minimum_distance_to_waters',
    'volume_overlap_with_protein',
    'volume_overlap_with_organic_cofactors',
    'volume_overlap_with_inorganic_cofactors',
    'volume_overlap_with_waters',
    'rmsd_≤_2å',
    'passes_valence_checks',
    'passes_kekulization',
    'inchi_crystal_valid',
    'inchi_docked_valid',
    'inchi_crystal',
    'inchi_docked',
    'inchi_overall',
    'inchi_version',
    'stereochemistry_preserved',
    'hydrogens',
    'net_charge',
    'protons',
    'stereo_sp3',
    'stereo_sp3_inverted',
    'stereo_type',
    'number_bonds',
    'shortest_bond_relative_length',
    'longest_bond_relative_length',
    'number_short_outlier_bonds',
    'number_long_outlier_bonds',
    'number_angles',
    'most_extreme_relative_angle',
    'number_outlier_angles',
    'number_noncov_pairs',
    'shortest_noncovalent_relative_distance',
    'number_clashes',
    'number_valid_bonds',
    'number_valid_angles',
    'number_valid_noncov_pairs',
    'number_aromatic_rings_checked',
    'number_aromatic_rings_pass',
    'aromatic_ring_maximum_distance_from_plane',
    'number_double_bonds_checked',
    'number_double_bonds_pass',
    'double_bond_maximum_distance_from_plane',
    'ensemble_avg_energy',
    'mol_pred_energy',
    'energy_ratio',
    'smallest_distance_protein',
    'num_pairwise_clashes_protein',
    'most_extreme_ligand_atom_id_protein',
    'most_extreme_protein_atom_id_protein',
    'most_extreme_ligand_element_protein',
    'most_extreme_protein_element_protein',
    'most_extreme_ligand_vdw_protein',
    'most_extreme_protein_vdw_protein',
    'most_extreme_sum_radii_protein',
    'most_extreme_distance_protein',
    'most_extreme_sum_radii_scaled_protein',
    'most_extreme_relative_distance_protein',
    'most_extreme_clash_protein',
    'smallest_distance_organic_cofactors',
    'not_too_far_away_organic_cofactors',
    'num_pairwise_clashes_organic_cofactors',
    'most_extreme_ligand_atom_id_organic_cofactors',
    'most_extreme_protein_atom_id_organic_cofactors',
    'most_extreme_ligand_element_organic_cofactors',
    'most_extreme_protein_element_organic_cofactors',
    'most_extreme_ligand_vdw_organic_cofactors',
    'most_extreme_protein_vdw_organic_cofactors',
    'most_extreme_sum_radii_organic_cofactors',
    'most_extreme_distance_organic_cofactors',
    'most_extreme_sum_radii_scaled_organic_cofactors',
    'most_extreme_relative_distance_organic_cofactors',
    'most_extreme_clash_organic_cofactors',
    'smallest_distance_inorganic_cofactors',
    'not_too_far_away_inorganic_cofactors',
    'num_pairwise_clashes_inorganic_cofactors',
    'most_extreme_ligand_atom_id_inorganic_cofactors',
    'most_extreme_protein_atom_id_inorganic_cofactors',
    'most_extreme_ligand_element_inorganic_cofactors',
    'most_extreme_protein_element_inorganic_cofactors',
    'most_extreme_ligand_vdw_inorganic_cofactors',
    'most_extreme_protein_vdw_inorganic_cofactors',
    'most_extreme_sum_radii_inorganic_cofactors',
    'most_extreme_distance_inorganic_cofactors',
    'most_extreme_sum_radii_scaled_inorganic_cofactors',
    'most_extreme_relative_distance_inorganic_cofactors',
    'most_extreme_clash_inorganic_cofactors',
    'smallest_distance_waters',
    'not_too_far_away_waters',
    'num_pairwise_clashes_waters',
    'most_extreme_ligand_atom_id_waters',
    'most_extreme_protein_atom_id_waters',
    'most_extreme_ligand_element_waters',
    'most_extreme_protein_element_waters',
    'most_extreme_ligand_vdw_waters',
    'most_extreme_protein_vdw_waters',
    'most_extreme_sum_radii_waters',
    'most_extreme_distance_waters',
    'most_extreme_sum_radii_scaled_waters',
    'most_extreme_relative_distance_waters',
    'most_extreme_clash_waters',
    'volume_overlap_protein',
    'volume_overlap_organic_cofactors',
    'volume_overlap_inorganic_cofactors',
    'volume_overlap_waters',
    'rmsd',
    'kabsch_rmsd',
    'centroid_distance']
    true_poses = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95_datamol_af2'
    pred_poses = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/results/cskde95_inference'
    rmsd_results = rmsd_only(true_poses, pred_poses)
    print(f'top1_rmsd_perc: {rmsd_results[0]}, median_top1: {rmsd_results[1]}, std_top1: {rmsd_results[2]}', f'\n{len(rmsd_results[3])}', ' proteins did not have poses')

if __name__ == "__main__":
    main()   