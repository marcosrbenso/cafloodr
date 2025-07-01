

#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
#   -----  DESCRIPTION  -----
#   Name : Forward-Mole 1D
#   Objetive: Correct DEM in urban areas according to a river shapefile in order to
#   improve hydrodynamic connectivity in water bodies
#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§



#alpha = 0.1  # This aids to control height differences considered to scan the channel width
#elements_ahead = 8  # This controls the number of rivers ahead considered to find a reference point
# --- Reading the required data --- #
#path2 = "C:/Projetos/13_Processamento_DEM/"
#dem_name = "dem_andadem"
#hydro_name = "strahler"


def FM1D(path2,dem_name,hydro_name,alpha = 0.1,elements_ahead = 8):
    # --- required libraries --- ##
    import os
    import time
    import rasterio
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import LineString

    import matplotlib
    matplotlib.use('TkAgg',force=True)
    import matplotlib.pyplot as plt

    start = time.time()

    #§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
    # Reading the original DEM
    with rasterio.open(path2+dem_name+'.tif') as src:
        dem_o = src.read(1)
        kwargs = src.meta
        bounds = src.bounds
        res = np.round(src.res[1])
        res2 = (res ** 2 + res ** 2) ** 0.5
        res_spacing = 0.001 * res
        dem_new = dem_o.copy()

    # The new raster name

    new_name = dem_name + '_1D_FM' + '.tif'
    


    # Load the shapefile
    gdf = gpd.read_file(path2+hydro_name+'.shp')
    line_bounds = gdf.total_bounds
    main = np.full((0, 7), np.nan, dtype=np.float32)
    lines = gdf.geometry
    df = pd.DataFrame()

    print('###### Updating the shapefile with row_id and strahler order ######')

    gdf = gdf[['FID','id_outlet',"STRAHLER","geometry"]]
    #gdf = gdf.rename(columns={gdf.columns[0]: 'id_outlet'})
    gdf = gdf.rename(columns={gdf.columns[2]: 'strahler'})
    gdf = gdf.rename(columns={gdf.columns[0]: 'row_id'})
    gdf["strahler"] = gdf["strahler"]
    gdf["row_id"] = gdf["row_id"]
    id_outlet = gdf['id_outlet']

    # Extracting the coordinates with their id element
    for i in range(len(lines)):
        if lines[i] is None:
            continue
        df_temp = pd.concat([pd.DataFrame(np.array(np.linspace(i,i,len(lines[i].coords[:])))) , pd.DataFrame(np.array(lines[i].coords[:]))],axis=1)
        df = pd.concat([df, df_temp], axis=0)
    df = np.array(df)


    print('###### Iterate over the shapefile features ######')
    for geom, strahler, row_id in zip(gdf.geometry, gdf['strahler'], gdf['row_id']):
        if geom and not geom.is_empty and geom.is_valid and geom.geom_type == 'LineString':
            # Convert line to points
            k = 1
            for i in range(len(geom.coords)-1):
                num_points = int(
                    ((geom.coords[i][0] - geom.coords[i+1][0]) ** 2 + (geom.coords[i][1] - geom.coords[i+1][1]) ** 2) ** 0.5 / res_spacing) + 1
                coordinates_array = np.linspace(geom.coords[i], geom.coords[i+1], num_points)
                coordinates_array = np.array([np.round((coordinates_array[:,0]-(bounds[0]-res/2))/res)*res+bounds[0]-res/2,
                                     np.round((coordinates_array[:,1]-(bounds[1]-res/2))/res)*res+bounds[1]-res/2]).T
                coordinates_array = coordinates_array[np.sort(np.unique(coordinates_array,axis=0, return_index=True)[1])]
                # Dropping the last row to avoid repeat the first and last point between segments
                coordinates_array = coordinates_array[:-1]
                # Populate the matrix with line points and additional data

                col = np.array((coordinates_array[:, 0] - bounds[0]) / res, dtype=int)
                row = np.array((bounds[3] - coordinates_array[:, 1]) / res, dtype=int)
                temp = np.vstack([col, row, np.repeat(strahler, len(col)), np.repeat(row_id, len(col)),
                                  coordinates_array[:, 0], coordinates_array[:, 1], range(k, k + len(col), 1)]).T
                main = np.vstack([main, temp])
                k = k + len(col)

    print('###### Overall forward mole loop ######')
    for stra in range(1, int(np.max(main[:,2]))+1):
        qq = 0
        if stra == np.max(main[:,2]):
            order = np.unique(main[np.argwhere(main[:, 2] == stra), 3]).astype(int)
            order = order[::-1]
        else:
            order = np.unique(main[np.argwhere(main[:, 2] == stra), 3]).astype(int)
        for j in order:
            qq += 1
            print(' ####################################################')
            print(' ##### Processing element from strahler order -' + str(stra) + '- #####')
            print(' ##### Processing element No. ' + str(j) + ' ### ' +
                  str(len(np.unique(main[np.argwhere(main[:, 2] == stra), 3]).astype(int)) - qq) + ' Elements left #####')
            # Segments to analyze
            temp_main = main[np.argwhere(main[:, 3] == j)[:, 0], :]

            if len(temp_main) == 0:  # if some line_id does not exist, the code will jump the line
                continue

            fl, lines_forward = j, 0
            # Creating a 3x3 empty windows
            matrix = np.full((3, 3, 3), np.nan, dtype=np.float32)
            t2 = np.zeros((8, 2))

            while lines_forward < elements_ahead:  # number of lines to see forward
                # Taking the last pair of coordinates (row and col)
                x_coor = temp_main[-1, 0]
                y_coor = temp_main[-1, 1]
                
                print(y_coor)

                # Window to identify neighboring coordinates
                t2[:, 0] = [x_coor-1, x_coor, x_coor+1, x_coor-1, x_coor+1, x_coor-1, x_coor, x_coor+1]
                t2[:, 1] = [y_coor-1, y_coor-1, y_coor-1, y_coor, y_coor, y_coor+1, y_coor+1, y_coor+1]
                # Finding the rows in main where required data is stored for our 3x3 window
                t3 = np.argwhere([(main[:, 0] == t2[value, 0]) & (main[:, 1] == t2[value, 1]) for value in range(len(t2))])[:, 1]
                # Filling the 3x3 window with rowd_id values
                matrix[0,(main[t3, 0] - min(t2[:, 0])).astype(int), (main[t3, 1] - min(t2[:, 1])).astype(int)] = main[t3, 3]
                # Filling the 3x3 window with k order line
                matrix[1,(main[t3, 0] - min(t2[:, 0])).astype(int), (main[t3, 1] - min(t2[:, 1])).astype(int)] = main[t3, 6]
                # Filling the 3x3 window with stralher line
                matrix[2,(main[t3, 0] - min(t2[:, 0])).astype(int), (main[t3, 1] - min(t2[:, 1])).astype(int)] = main[t3, 2]
                # Dropping not allowed neighbors
                # Dropping Nan values
                # Transforming t2 values of rows and columns for a window of 3x3 size
                t4 = t2.copy()
                t4[:,0], t4[:,1] = t4[:,0]-min(t4[:,0]), t4[:,1]-min(t4[:,1])
                t4 = t4.astype(int)
                condition = np.isnan(matrix[0, t4[:, 0], t4[:, 1]])

                # Windows of the row_id values
                condition2 = np.array([value in main[:,3] for value in matrix[0, t4[:, 0], t4[:, 1]]])

                # Dropping values that are lower than the last value K and belongs to the same row_id
                condition2[np.argwhere(condition2==True)[((matrix[1,t4[condition2][:,0],t4[condition2][:,1]]<temp_main[-1,6])&
                                                         (matrix[0,t4[condition2][:,0],t4[condition2][:,1]]==temp_main[-1,3]))]] = False

                # Dropping values that are lower than the last value K and with lower stralher order
                t4 = np.where(condition.reshape(-1, 1), np.nan, t4)
                t4 = np.where(condition2.reshape(-1, 1), t4, np.nan)
                t4 = t4[~np.isnan(t4).any(axis=1)].astype(int)

                # sorting t2 indexes
                t4 = t4[np.lexsort((matrix[1,t4[:,0],t4[:,1]], matrix[0,t4[:,0],t4[:,1]]))]

                # returning original values of t4 for row and columns as t2
                t4[:,0], t4[:,1] = t4[:,0]+min(t2[:,0]), t4[:,1]+min(t2[:,1])

                # if t4 is empty, means that could be the outlet because there is no more parts downstream of it
                if len(t4)==0:
                    break
                # extracting all info from main for the filtered cells
                temp2 = main[np.argwhere([(main[:,0]==t4[value,0])&(main[:,1]==t4[value,1]) for value in range(len(t4))])[:,1],:]

                if len(temp2) == 0: # this check if temp 2 is empty, this happend when we are working with the outlet line (last one)
                    break
                if np.any(temp2[:, 3] != temp2[0, 3]):
                    # if there is row_id different, we sort according x[2] and x[6] which are stralher number and k order, respectively.
                    if np.all(temp2[:,2]==temp2[0,2]):
                        # Finding the group with the largest K order, this means that is the end of one channel, we looking for
                        # the beginning of the downstream channel
                        df = pd.DataFrame(temp2)
                        # Ordening the pandas group by the mean mentioned before.
                        temp2 = pd.concat([df[df[3] == df.groupby(3)[6].mean().idxmax()],df[df[3] != df.groupby(3)[6].mean().idxmax()]]).to_numpy()
                    else:
                        df = pd.DataFrame(temp2)
                        grouped_df = df.groupby(3)[6].mean()
                        sorted_group_indices = np.argsort(df.groupby(3)[6].mean().values)[::-1]
                        # Sort the DataFrame based on the sorted group indices
                        temp2 = df.merge(grouped_df, left_on=3, right_index=True).sort_values(by='6_y', ascending=False).to_numpy()[:, :-1]
                
                max_ref_attempts = 3
                
                
                
                
                if temp2[-1, 2] >= temp_main[-1, 2]:
                    # gather downstream segment candidates
                    t = main[ main[:, 3] == temp2[-1, 3] ]
                    if len(t) <= 2:
                        temp_main = np.vstack([temp_main, t])
                    else:
                        # try at most max_ref_attempts times to find a good matching reference
                        for attempt in range(max_ref_attempts):
                            # pick the "best" candidate in temp2
                            idxs = np.argwhere(
                                (temp2[:,2] == temp2[:,2].max()) &
                                (temp2[:,6] == temp2[:,6].min())
                            ).flatten()
                
                            if idxs.size == 0:
                                print("Warning: no matching reference coordinate in temp2.")
                                break
                
                            # reference coordinate
                            rc = temp2[idxs[0], :2]
                            # find it in t
                            t_idx = np.argwhere((t[:,0] == rc[0]) & (t[:,1] == rc[1])).flatten()
                            if t_idx.size == 0:
                                print(f"Warning: ref {rc} not in t; attempt {attempt+1}")
                                continue   # try again
                            # success!
                            t = t[t_idx[0]:, :]
                            temp_main = np.vstack([temp_main, t])
                            break
                        else:
                            # if we exhausted all attempts without break:
                            print(f"Failed to match reference after {max_ref_attempts} attempts—skipping element.")
                            flag_out_element = 1

                else:
                    temp_main = np.vstack([temp_main, temp2])
                lines_forward += 1
                fl = temp_main[-1,6].astype(int)

            # Calculate distances between points according to the coordinates
            temp = temp_main
            temp = np.hstack([temp, np.zeros((temp.shape[0], 1))])
            for i in range(1, len(temp)):
                temp[i, -1] = np.sqrt((temp[i, 4]-temp[i-1, 4])**2 + (temp[i, 5]-temp[i-1, 5])**2)

            flag_r = 0 # to indicate that right comes first
            flag_l = 0
            flag_2D_off = 0
            flag_2nd_loop_out = 0
            flag_2D_repeat = 0
            flag_backup = 0
            flag_out_element = 0
            flag_first = 0
            flag_once = 0
            k_backup = 0
            flag_empty = 0
            flag_comment = 0

            t = -1
            while t < 10000:
                flag_2nd_loop_out = 0

                if flag_out_element == 1:
                    break

                if t == -1:  # the first iteration for the center line
                    temp_w = temp.copy()
                    if flag_empty == 1:
                        flag_2D_repeat = 0
                else:
                    if flag_2D_repeat == 1 and flag_comment == 0:
                        print('...fixing the right side')
                        flag_comment = 1

                    if t > np.size(width_right, 0)-1 and flag_r == 0:
                        print('...fixing the left side')
                        flag_r = 1
                        if flag_once == 0:
                            flag_first = 1
                            flag_once = 1  # only once per segment
                        t = 0  # now is left turn

                    if flag_r == 0:
                        temp_w = temp.copy()
                        temp_w[:, 0:2] = temp_w[:, 0:2] + width_right[t, :]
                        if flag_once == 1:
                            flag_empty = 1
                        t = t + 1
                    else:
                        temp_w = temp.copy()
                        if t > np.size(width_left, 0)-1:
                            temp_w = temp.copy()  #  here we return to the central line
                            flag_backup = 1
                            flag_r = 0
                            flag_2D_repeat = 0  # to activate the width scan
                            if flag_first == 1:
                                flag_first = 0
                            t = -1
                        else:
                            temp_w[:, 0:2] = temp_w[:, 0:2] + width_left[t, :]
                            t = t + 1

                # Main Forward Mole script
                i = j

                if flag_2D_repeat == 1 and flag_first == 1:
                    k = 0
                elif flag_2D_repeat == 0 and flag_empty == 1:
                    k = k_0
                    flag_empty = 0
                elif flag_2D_repeat == 1 and flag_empty == 1 and flag_once == 1:
                    k = k_0
                    flag_empty = 1
                else:
                    k = 0

                if flag_backup == 1:
                    i = i_backup
                    k = k_backup
                    flag_backup = 0
                    
                        
                max_pit_attempts = 3  # Maximum attempts to fix the same pit
                pit_attempts = 0
                
                max_k = len(temp_w)
            

                while i == j and flag_2nd_loop_out == 0:
                    print(k)
                    # Check for infinite loop condition first
                    if pit_attempts > max_pit_attempts:
                        print(f"Max pit correction attempts reached at k={k}. Moving forward.")
                        k += 1
                        if k <= max_k:
                            flag_out_element = 1
                        pit_attempts = 0
                        continue
                    
                    # at the top of your inner pit‐detection loop:
                    
                    elif k >= max_k - 1:
                        flag_out_element = 1
                        
                        break
                    
                    elif dem_o[int(temp_w[k,1]), int(temp_w[k,0])] < dem_o[int(temp_w[k+1, 1]), int(temp_w[k+1, 0])]:

                        print(f'-------- Pit found at {k*res} meters')
                        pit_attempts += 1  # Increment attempt counter
                        

                        # Check for valid solutions
                        valid_solutions = np.argwhere((dem_o[int(temp_w[k, 1]), int(temp_w[k, 0])] >= dem_o[temp_w[:, 1].astype(int), temp_w[:, 0].astype(int)]) & 
                                          (dem_o[temp_w[:, 1].astype(int), temp_w[:, 0].astype(int)] < dem_o[int(temp_w[k,1]), int(temp_w[k,0])]))
                        

                        
                        if len(valid_solutions) == 0:
                            flag_out_element = 1
                            
                            break

                        

                        #######################################################################################################
                        # identify the cells that belongs to the error and solution by the FM logic
                        refs = np.argwhere((dem_o[int(temp_w[k, 1]), int(temp_w[k, 0])] > dem_o[temp_w[:, 1].astype(int), temp_w[:, 0].astype(int)]) &
                                                       (dem_o[temp_w[:, 1].astype(int),temp_w[:, 0].astype(int)] < dem_o[int(temp_w[k, 1]),int(temp_w[k, 0])]))
                        
                        
                        if len(refs)<=1:
                            k = k + 1
                            continue
                        
                        
                        elif len(np.argwhere(k < refs)) == 0:
                            k = k + 1
                            continue
                        refs = np.arange(k, refs[np.argwhere(k<refs)[0][0].astype(int)], 1)
                        
                        # and temp_w[refs[-1], 3] != i
                        if t != -1:
                            if refs[-1] >= refs_backup * 1:
                                # if j ==5 and flag_r == 1 and t == 8 and k > 700:
                                #     print('t = ' + str(t) + 'and k = ' + str(k))
                                if len(np.argwhere(refs < refs_backup)) == 0:  # no solution
                                    flag_2nd_loop_out = 1
                                    continue
                                refs = refs[:np.argwhere(refs == refs_backup)[-1][0].astype(int)+1]
                                # flag_2nd_loop_out = 1
                                if (t > np.size(width_left, 0)-1) and flag_r == 1:
                                    flag_r = 0
                                    flag_backup = 1
                                    t = -1  #  to get back to the central line

                                dem_o[int(temp_w[int(refs[-1]), 1]), int(temp_w[int(refs[-1]), 0])] = dem_backup
                                


                        if len(refs) == 0:
                            k = k + 1
                            if k > refs_backup and t != -1:
                                flag_2nd_loop_out = 1
                            continue

                        # for the identified cells, we apply linear interpolation.
                        
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        print("Interpolated")
                        dem_o[temp_w[refs[1:-1],1].astype(int),temp_w[refs[1:-1],0].astype(int)] = \
                            ((dem_o[int(temp_w[refs[0],1]),int(temp_w[refs[0],0])] - dem_o[int(temp_w[refs[-1], 1]),int(temp_w[refs[-1], 0])]) *
                             (np.cumsum(temp_w[refs[1:-1], 7]) - np.sum(temp_w[refs[1:], 7]))) / (
                                           temp_w[refs[0], 7] - np.sum(temp_w[refs[:], 7])) + dem_o[int(temp_w[refs[-1], 1]), int(temp_w[refs[-1], 0])]

                        #######################################################################################################


                        ######################################################################################################

                        if flag_2D_off == 1:

                            k = refs[-1]
                            k_0 = refs[0]
                            # i = temp_w[k, 3]
                            flag_2nd_loop_out = 1   # to get out from the center line and work the river width
                            flag_2D_off = 0  # to avoid lose k and i backup in the next iteration when river width is fixed
                            flag_comment = 0
                            i_backup = i
                            k_backup = k
                            dem_backup = dem_o[int(temp_w[k, 1]),int(temp_w[k, 0])]
                            refs_backup = refs[-1]
                            t = 0  # to get out from the center line
                            if flag_once == 0:
                                flag_first = 1
                        else:
                            if flag_2D_repeat == 1 and refs[-1] >= refs_backup:
                                flag_2nd_loop_out = 1
                                flag_empty = 1
                                if flag_r == 1 and flag_l == 1:
                                    t = -1  # to return t the central line
                                    flag_empty = 1
                                    flag_r = 0
                            else:
                                k = refs[-1]
                                i = temp_w[k, 3]
                    else:
                        if flag_2D_repeat == 1 and k >= refs_backup:
                            flag_2nd_loop_out = 1
                        else:
                            k += 1
                            i = temp_w[k, 3]

                    if flag_r == 1 and t > np.size(width_left, 0)-1:
                        i = temp_w[k_backup, 3]
                        flag_first = 0

                    if i != j:
                        flag_out_element = 1
                        
                        # After processing, FORCE k increment if still in same location
                        if not flag_2D_repeat:  # Only increment if not in 2D mode
                            k += 1
                            pit_attempts = 0
                            

                    else:
                        # No pit found, move forward
                        k += 1
                        pit_attempts = 0  # Reset counter                    

    # saving the new raster
    with rasterio.open(path2+"/"+new_name, 'w', **kwargs) as dst:
        dst.write_band(1, dem_o.astype(rasterio.float32))
        
    end = time.time()
    print('Modified DEM successfully exported')
    print('Time elapsed: ' + str((end-start)/60) + ' minutes')