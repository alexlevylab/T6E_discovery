from pathlib import Path
import json
import glob
import re
import numpy as np
import pandas as pd


def process_files(globpath):
    '''
    Parameters
    ----------
    globpath : str
        Path with wildcard of json files.

    Returns
    -------
    dictionary of files and paths with parsed names
    '''
    all_files = glob.glob(globpath)
    # print(all_files)
    
    xs = []
    ys = []
    filenames = {}
    
    for file in all_files:
        scores = json.loads(Path(file).read_text())
        x = re.search("(\d+)_and", file)[1]
        y = re.search("and_(\d+)_", file)[1]
        filenames[x + "__"+ y]  = file
        
        xs.append(x)
        ys.append(y)
        
    xs = list(set(xs))
    ys = list(set(ys))
    
    xs.sort()
    ys.sort()
    
    return filenames, xs, ys

    
def define_good_tracts(plddt):
    '''
    Parameters
    ----------
    plddt : DataFrame
        dataframe with plddt values
    
    Returns
    -------
    good_tracts : DataFrame
        DataFrame with info about tracts

    '''
    rolling_score_plddt = plddt[0].rolling(window = 30).mean().to_frame()
    rolling_score_plddt[1] = rolling_score_plddt[0] > 70
    rolling_score_plddt.reset_index(inplace = True)
        
    good_tracts = rolling_score_plddt.loc[rolling_score_plddt[1] == True]        
    good_tracts['shifted_index_f'] = good_tracts['index'].shift(1)
    good_tracts['diff_f']  =  good_tracts['index'] - good_tracts['shifted_index_f']
    good_tracts['spaced_f'] = good_tracts['diff_f'] == 1
    good_tracts['gap_size'] = good_tracts['spaced_f']
    good_tracts['tract_ID'] = good_tracts.loc[good_tracts['gap_size']==False].groupby('gap_size').cumcount() 
    good_tracts['tract_ID'] = good_tracts['tract_ID'].fillna(method='ffill')
    
    return good_tracts
    

def calculate_rectangle_reversed_normalized_mean(pae, start1, end1, start2, end2):
    '''
    Takes two rectangles in pae matrix, normalize score between 0 and 1 (1 being best confidence)
    These rectangles represent the -,- and +,- quadrants, where we expect to see signs of intermolecular PAE scores
    
    Parameters
    ----------
    pae : DataFrame
        pae matrix of scores from 0-30

    start :
        start coordinate of interest of square
        
    end : 
        end coordinate of interest of square 

    Returns
    -------
    mean and std of flattened score matrix in intermolecular region
    '''
    reversed_normalized_pae = (30 - pae)/30
    rec1 = reversed_normalized_pae.iloc[start1:end1, start2:end2]
    rec2 = reversed_normalized_pae.iloc[start2:end2, start1:end1]

    rec1_flat = rec1.to_numpy().flatten()
    rec2_flat = rec2.to_numpy().flatten()
 
    allscores_flat = np.concatenate((np.array(rec1_flat), np.array(rec2_flat)))

    return allscores_flat.mean(), allscores_flat.std()

def scoring(x,y, filenames, df):
    '''
    Parameters
    ----------
    x : str
        processed name of toxin that is a key in dictionary of files (toxin name)

    y : str
        processed name of toxin that is a key in dictionary of files (imm name)
        
    filenames : dict
        dictionary with filenames. x and y are keys to this dict
        
    df : DataFrame
        dataframe with size info

    Returns
    -------
    scoring
    '''
    
    out = []

    # load json, save in df
    scores = json.loads(Path(filenames[x + "__" + y]).read_text())
    plddt = pd.DataFrame(scores['plddt'])
    pae = pd.DataFrame(scores['pae'])
    
    # calculate labels and lengths
    label = x + " and " + y
    length_x = int(df.loc[df['Gene ID'] == x]['Amino Acid Sequence Length (aa)'].iloc[0])
    length_y = int(df.loc[df['Gene ID'] == y]['Amino Acid Sequence Length (aa)'].iloc[0])
    length_together = length_x + length_y

    #are they expected?
    digitsX = [int(s) for s in x.split() if s.isdigit()]
    digitsY = [int(s) for s in y.split() if s.isdigit()]

    expected = abs(digitsX[0] - digitsY[0]) <= 1
    print(expected)
    # call function to find good tracts
    good_tracts_x = define_good_tracts(plddt.iloc[range(0, length_x)])
    good_tracts_y = define_good_tracts(plddt.iloc[range(length_x, length_together)])
    
    
    if (not good_tracts_x.empty) and (not good_tracts_y.empty):
    
        # find start coordinates of good tracts                        
        starts_x = good_tracts_x.groupby('tract_ID')['index'].min().to_frame().rename(columns = {'index' : 'start_coord'})
        starts_y = good_tracts_y.groupby('tract_ID')['index'].min().to_frame().rename(columns = {'index' : 'start_coord'})
        
        # calculate end coordinates of good tracts
        starts_x['end_coord'] = starts_x['start_coord'].shift(-1).fillna(good_tracts_x['index'].iloc[-1]).astype(int)
        starts_y['end_coord'] = starts_y['start_coord'].shift(-1).fillna(good_tracts_y['index'].iloc[-1]).astype(int)
        
        for i in starts_x.index.unique():
            
            start_coord_x_tract = int(starts_x.iloc[int(i)]['start_coord']) 
            end_coord_x_tract = int(starts_x.iloc[int(i)]['end_coord'])
            
            for j in starts_y.index.unique():
    
                start_coord_y_tract = int(starts_y.iloc[int(j)]['start_coord'])
                end_coord_y_tract = int(starts_y.iloc[int(j)]['end_coord'])
                
                
                
                intermolecular_score = calculate_rectangle_reversed_normalized_mean(pae, start_coord_x_tract, 
                                                                                     end_coord_x_tract, 
                                                                                     start_coord_y_tract, 
                                                                                     end_coord_y_tract)
                
                out.append([label, 
                            i, 
                            j, 
                            start_coord_x_tract, 
                            end_coord_x_tract, 
                            start_coord_y_tract, 
                            end_coord_y_tract,
                            intermolecular_score[0],
                            intermolecular_score[1],
                            expected,
                            digitsX[0],
                            digitsY[0]])
                
    else:
        out.append([label, 
                    None, 
                    None, 
                    None, 
                    None, 
                    None, 
                    None,
                    None,
                    None,
                    expected,
                    digitsX[0],
                    digitsY[0]])
                
    outframe = pd.DataFrame(out, columns = ['protein_names', 
                                            'tract ID top left', 
                                            'tract ID bottom right', 
                                            'start_coord_top_left_tract', 
                                            'end_coord_top_left_tract', 
                                            'start_coord_bottom_right_tract', 
                                            'end_coord_bottom_right_tract',
                                            'mean_intermolecular_PAE',
                                            'std_dev_intermolecular_PAE',
                                            'expected',
                                            'x_name',
                                            'y_name'])
                                
    return outframe
            
    
########################

# Import size data
df = pd.read_csv("sizes_and_metadata3.tsv", sep = '\t', dtype = 'str')

globpath = "./*unrelaxed_rank_1_model_*_scores.json"
filenames, xs, ys = process_files(globpath)

final = []

for file in filenames:
    x = re.search("(\d+)__", file)[1]
    y = re.search("__(\d+)", file)[1]       
    final.append(scoring(x, y, filenames, df))
        
final_df = pd.concat(final)
final_df.reset_index(drop = True, inplace = True)
final_df.fillna(value = 0, inplace = True)

idx = final_df.groupby("protein_names")['mean_intermolecular_PAE'].idxmax()

filtered_df = final_df.loc[idx]
filtered_df.sort_values(['x_name', 'expected'], inplace = True)

print(filtered_df)

filtered_df.to_csv("scores_maxima.csv", index = False)
