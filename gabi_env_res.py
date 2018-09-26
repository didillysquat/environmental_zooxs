import pandas as pd
import sys
import pickle
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import matplotlib
from collections import defaultdict
import sys
import statistics
import itertools
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import numpy as np
from multiprocessing import Queue, Process, Manager


def create_diverstiy_figs():
    ''' This will be the code to create the diversity figure for Gabi's paper.

        We have used symportal to do the QC work on the samples.
        As such we have tab delimited outputs that we can read into pandas dataframes
        We also have the info list from Gabi that contains what all of the samples are
        Esentially what we want to do is quite simple, however, we will need to be careful with the colouring
        As there are going to be a shed ton of sequences we will probably want to order the sequences found some how
        and colour those, then maybe colour sequences found at below 5% to grey. Something like that.
        I have a load of code from previous work when we were trying to look for actual SymPortal ITS2 type profiles
        within the env samples but we are no longer doing that. Still some of the code below might be useful.
        Actually I'm going to drop it into a junk code def and start from scratch here.
    '''

    try:
        sp_output_df = pickle.load(open('sp_output_df.pickle', 'rb'))
        QC_info_df= pickle.load(open('QC_info_df.pickle', 'rb'))
        info_df = pickle.load(open('info_df.pickle', 'rb'))
    except:
        # read in the SymPortal output
        sp_output_df = pd.read_csv('131_142.DIVs.absolute.txt', sep='\t', lineterminator='\n')

        # The SP output contains the QC info columns between the DIVs and the no_name ITS2 columns.
        # lets put the QC info columns into a seperate df.
        QC_info_df = sp_output_df[['Samples','raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique']]

        # now lets drop the QC columns from the SP output df and also drop the clade summation columns
        # we will be left with just clumns for each one of the sequences found in the samples
        sp_output_df.drop(columns=['noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                                'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                                'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique'
                                   ]
                          , inplace=True)

        # read in the info file
        info_df = pd.read_csv('info_300718.csv')

        # drop the rows that aren't from the 22.08.2016 data
        info_df = info_df[info_df['date collected'] == '22.08.2016']

        # Now we need to link the SP output to the sample names in the excel. Annoyingly they are formatted
        # slightly differently so we can't make a direct comparison.
        # easiest way to link them is to see if the first part of the SP name is the same as the first part
        # of the 'sequence file' in the meta info
        # when doing this we can also drop the SP info for those samples that won't be used i.e. those that
        # aren't now in the info_df

        # firstly rename the colums so that they are 'sample_name' in all of the dfs
        QC_info_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        sp_output_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        info_df.rename(index=str, columns={'Sample Name': 'sample_name'}, inplace=True)

        indices_to_drop = []
        for sp_index in sp_output_df.index.values.tolist():
            # keep track of whether the sp_index was found in the info table
            # if it wasn't then it should be dropped
            found = False
            for info_index in info_df.index.values.tolist():
                if sp_output_df.loc[sp_index, 'sample_name'].split('_')[0] == info_df.loc[info_index, 'Sequence file'].split('_')[0]:
                    found = True
                    # then these are a related set of rows and we should make the sample_names the same
                    sp_output_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']
                    QC_info_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']


            if not found:
                indices_to_drop.append(sp_index)

        # drop the rows from the SP output tables that aren't going to be used
        sp_output_df.drop(inplace=True, index=indices_to_drop)
        QC_info_df.drop(inplace=True, index=indices_to_drop)

        # let's sort out the 'environ type' column in the info_df
        # currently it is a bit of a mess
        for index in info_df.index.values.tolist():
            if 'coral' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'coral'
            elif 'seawater' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sea_water'
            elif 'mucus' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'mucus'
            elif 'SA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_close'
            elif 'SB' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_far'
            elif 'TA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'turf'

        # now clean up the df indices
        info_df.index = range(len(info_df))
        sp_output_df.index = range(len(sp_output_df))
        QC_info_df.index = range(len(QC_info_df))

        # make the sample_name column the index for each of the datasets
        info_df.set_index('sample_name', inplace=True)
        sp_output_df.set_index('sample_name', inplace=True)
        QC_info_df.set_index('sample_name', inplace=True)

        # pickle the out put and put a check in place to see if we need to do the above
        pickle.dump(sp_output_df, open('sp_output_df.pickle', 'wb'))
        pickle.dump(QC_info_df, open('QC_info_df.pickle', 'wb'))
        pickle.dump(info_df, open('info_df.pickle', 'wb'))



    # so we want to plot the ITS2 sequence diversity in each of the samples as bar charts
    # We are going to have a huge diversity of sequences to deal with so I think something along the lines
    # of plotting the top n most abundant sequences. The term 'most abundant' should be considered carefully here
    # I think it will be best if we work on a sample by sample basis. i.e. we pick the n sequences that have the
    # highest representation in any one sample. So for example what we are not doing is seeing how many times
    # C3 was sequenced across all of the samples, and finding that it is a lot and therefore plotting it.
    # we are looking in each of the samples and seeing the highest proportion it is found at in any one sample.
    # This way we should have the best chance of having a coloured representation for each sample's most abundant
    # sequence.

    # to start lets go sample by sample and see what the highest prop for each seq is.

    # dict to hold info on which sample and what the proportion is for each sequence
    # key = sequence name, value = tup ( sample name, relative abundance)
    try:
        seq_rel_abund_calculator_dict = pickle.load(open('seq_rel_abund_calculator_dict.pickle', 'rb'))
    except:
        seq_rel_abund_calculator_dict = {}
        for sample_index in sp_output_df.index.values.tolist():
            sys.stdout.write('\nGeting rel seq abundances from {}\n'.format(sample_index))
            # temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum(axis='index'))
            temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum())
            for seq_name in temp_prop_array.keys():
                sys.stdout.write('\rseq: {}'.format(seq_name))
                val = temp_prop_array[seq_name]
                if val != 0:  # if the sequences was found in the sample
                    # If the sequence is already in the dict
                    if seq_name in seq_rel_abund_calculator_dict.keys():
                        # check to seee if the rel abundance is larger than the one already logged
                        if val > seq_rel_abund_calculator_dict[seq_name][1]:
                            seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
                    # if we haven't logged for this sequence yet, then add this as the first log
                    else:
                        seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
        pickle.dump(seq_rel_abund_calculator_dict, open('seq_rel_abund_calculator_dict.pickle', 'wb'))



    # here we have a dict that contains the largest rel_abundances per sample for each of the seqs
    # now we can sort this to look at the top ? n sequences to start with (I'm not sure how the colouring will
    # look like so lets just start with 30 and see where we get to)
    sorted_list = sorted(seq_rel_abund_calculator_dict.items(), key = lambda x: x[1][1], reverse=True)
    most_abund_seq_names = [tup[0] for tup in sorted_list]

    # from the above sorted list we can then plot these sequences with colour and all others with grey scale
    # lets make a coulour dictionary for the most common types
    colour_list = get_colour_list()
    colour_dict = {}
    num_coloured_seqs = 30

    # we will also need a grey palette for those sequences that are not the ones being annotated
    grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    # make the dict
    for i in range(len(most_abund_seq_names)):
        if i < num_coloured_seqs:
            colour_dict[most_abund_seq_names[i]] = colour_list[i]
        else:
            colour_dict[most_abund_seq_names[i]] = grey_palette[i % 6]
    # add the 'low' key and assign it to grey for later on
    # the low category will be created later on to hold the grouped abundances of sequences in samples
    # below a certain rel abund cutoff
    colour_dict['low'] = '#D0CFD4'

    # for plotting we will also need to collect the 'rectangle' objects from the plotting process to use to make
    # the lengend which we will display on the plot right at the end
    # this is going to be a little tricky to collect as not ever sample type group that we are plotting
    # is going to have all of the top n sequences. So, we will have to pick up rectangles when they come up in
    # the plotting.
    # to keep track of which sequences we need we have:
    top_n_seqs_to_get = most_abund_seq_names[:num_coloured_seqs]
    # to keep track of which sequences we have already collected we have:
    top_n_seqs_done = []
    #finally to collect the rectangles and store them in order despite collecting them out of order we have:
    legend_rectangle_holder = [[] for i in range(num_coloured_seqs)]


    # set up the plotting environment
    # one plot for coral, mucus, seawater, sediment, turf_algae
    f, axarr = plt.subplots(6, 1)
    # counter to reference which set of axes we are plotting on
    axarr_index = 0

    # now we need to get the actual plotting information for each of the samples.
    # we can do this sample type by sample type
    # we will create a local dataframe that will be a sub set of the main sp output dataframe but will
    # only contain the samples of the given sample type. It will also eventually only contain the sequence
    # information for the sample type in question. Normally I would make a 2D list to hold the plotting
    # information but in this case having all of the informatin in a dataframe is working really well and
    # this is how I will do it in future.

    # go environment type by environment type
    for env_type in ['coral', 'mucus', 'sea_water', 'sed_close', 'sed_far', 'turf']:
    # for env_type in ['sed_far',  'turf']:
        sys.stdout.write('\nGenerating plotting info for {} samples\n'.format(env_type))
        # currently we have something like 4000 sequences to plot which is too may
        # I think it will be much easier if we group the sequences that are found below a certain
        # threshold. I think the best way to do this is to create a slice of the main df that
        # contain the information for the samples of the env_type only

        # get subset of the main dfs that contain only the coral samples
        env_info_df = info_df[info_df['environ type'] == env_type]
        env_sp_output_df = sp_output_df.loc[env_info_df.index.values.tolist()]

        # append a 'low' columns to the env_sp_ouput_df and populate with 0s
        env_sp_output_df['low'] = 0
        # now make a proportions version of the df, rather than absolute counts
        env_sp_output_df_prop = env_sp_output_df[:].div(env_sp_output_df[:].sum(axis=1), axis=0)

        # now as we work our way through we can sum up the low sequences into this column
        # we can then check for 0 columns and drop these.

        # get a list of the sequences found in the collection of samples of the given type and order
        # them according to summed rel_abundance acorss all samples. This should be the order in which
        # the samples are plotted
        # at the same time we can get the info we need for plotting

        summed_seq_rel_abund_across_smpls_dict = {seq: 0 for seq in list(sp_output_df)}

        # we are also going to need to put the samples in an order that makes sense.
        # Ideally I would like to take the time to do a proper travelling salesman analysis
        # and I was thinking of putting in an ant colony system to do this but...
        # I'm not sure we've got time for that. So let's do something a little more simple which
        # should still be effective and look groovy
        # Lets sort according to the majority sequence. I.e. lets put the samples into groups that
        # are defined by what their majorty sequence is, then lets plot in order of those groups.
        # within the groups we will plot in the order of the most abundant rel abund first.
        # to get this information we will have a dict to collect it. The key of the dict will
        # be the majority sequence with the value being a list that contains tuples. One tuple
        # for each sample in the list which will contain the sample_name and the re abund of the maj seq

        # Dict for holding the sample_sorting info
        sample_sorting_info_dict = defaultdict(list)

        for sp_index in env_sp_output_df_prop.index.values.tolist():
            sys.stdout.write('\r{}'.format(sp_index))
            # we need to get the name of the most abundant sequence and its rel abund for each sample
            # for the sample_sorting_info_dict
            most_abund_seq_name = env_sp_output_df_prop.loc[sp_index].idxmax(axis='index')
            rel_abund_of_most_abund_seq = env_sp_output_df_prop.loc[sp_index, most_abund_seq_name]
            sample_sorting_info_dict[most_abund_seq_name].append((sp_index, rel_abund_of_most_abund_seq))
            # Put its seq rel abundances ino the summed_seq_rel_abund_across... dict

            for non_zero_seq in env_sp_output_df_prop.loc[sp_index][env_sp_output_df_prop.loc[sp_index] > 0].index:  # not including the final 'low' columns (this will be zero)
                val = env_sp_output_df_prop.loc[sp_index, non_zero_seq]
                # be sure to count the value of this cell and using it in judging which are the most
                # abundant sequences before we check whether to relegate it to the 'low' column
                summed_seq_rel_abund_across_smpls_dict[non_zero_seq] += val
                if val < 0.005:
                    env_sp_output_df_prop.loc[sp_index, 'low'] += val
                    env_sp_output_df_prop.loc[sp_index, non_zero_seq] = 0


        # here we can get a sorted sample list using the sample_sorting_info_dict
        sorted_sample_list = []
        # we want to work through the sample_sorting_info_dict by the longest lists first
        sorted_keys_for_sample_sort = \
            [tup[0] for tup in sorted(sample_sorting_info_dict.items(), key=lambda x: len(x[1]), reverse=True)]
        # for each of teh maj seq groups
        for sorted_key in sorted_keys_for_sample_sort:
            # now within each of these lists we want to order according to the rel_abundance of the sequences
            sorted_list_of_samples_of_group = [tup[0] for tup in sorted(sample_sorting_info_dict[sorted_key], key=lambda x: x[1], reverse=True)]
            sorted_sample_list.extend(sorted_list_of_samples_of_group)

        # now we should re-order the df so that it is in the sample order of sorted_sample_list
        env_sp_output_df_prop = env_sp_output_df_prop.reindex(sorted_sample_list)

        # here we have a dict that contains the abundances of the seqs for the coral samples
        # we will plot the coral samples' seqs in the order of the sequences in this dict
        sorted_list_of_env_specific_seqs_tup \
            = sorted(summed_seq_rel_abund_across_smpls_dict.items(), key=lambda x: x[1], reverse=True)
        sorted_list_of_env_specific_seqs = [tup[0] for tup in sorted_list_of_env_specific_seqs_tup]

        # Now we check for zero only columns and drop them
        # we also need to remove these columns from the sorted
        non_zero_cols = list(env_sp_output_df_prop.loc[:, (env_sp_output_df_prop != 0).any(axis=0)])
        sorted_list_of_env_specific_seqs = [seq for seq in sorted_list_of_env_specific_seqs if seq in non_zero_cols]

        # add the 'low' column
        sorted_list_of_env_specific_seqs.append('low')

        # now drop the cols
        env_sp_output_df_prop = env_sp_output_df_prop[non_zero_cols]

        # we also have the plotting list which contains the info we will be plotting.

        # the plotting_list is currently organised in a different order to that of the sorted_list_of_env...
        # we need to change this order
        env_sp_output_df_prop = env_sp_output_df_prop[sorted_list_of_env_specific_seqs]

        # we now need to transpose this
        env_sp_output_df_prop = env_sp_output_df_prop.transpose()

        # here we finally have the plotting lists in the order of the sorted_list_of_env
        bottom = [0 for smp in list(env_sp_output_df_prop)]
        bar_list = []
        # for each sample
        plotting_indices = range(len(list(env_sp_output_df_prop)))

        # for each sequence
        list_of_seqs = env_sp_output_df_prop.index.values.tolist()
        for i in range(len(list_of_seqs)):

            sys.stdout.write('\r{}/{}'.format(i, len(list(env_sp_output_df_prop.iloc[:, 0]))))
            bar_list.append(axarr[axarr_index].bar(plotting_indices, list(env_sp_output_df_prop.iloc[i]), 1, bottom,
                                                   color=colour_dict[sorted_list_of_env_specific_seqs[i]]))
            bottom = [L + M for L, M in zip(bottom, list(env_sp_output_df_prop.iloc[i]))]

            # check to see if the seq we are plotting is still in the top_n_seqs list
            # if it is in this list then we still need to grab a rectangle for plotting
            # the legend. Once we have grabbed the rectangle then we should remove the seq from the top_n_seqs list
            seq_name = list_of_seqs[i]

            if seq_name in top_n_seqs_to_get and seq_name not in top_n_seqs_done:
                # then this is a seq that we still need to get a rectangle for the legend
                legend_rectangle_holder[top_n_seqs_to_get.index(seq_name)].append(bar_list[-1][0])
                top_n_seqs_done.append(seq_name)

        # https://stackoverflow.com/questions/12998430/remove-xticks-in-a-matplotlib-plot
        axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        # https://stackoverflow.com/questions/15858192/how-to-set-xlim-and-ylim-for-a-subplot-in-matplotlib
        axarr[axarr_index].set_xlim((-0.5, (len(list(env_sp_output_df_prop)) - 0.5)))
        # https://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib
        axarr[axarr_index].spines['right'].set_visible(False)
        axarr[axarr_index].spines['top'].set_visible(False)

        # To make the legend for this mo fo is going to be a little tricky. The legend argument basically takes
        # two lists. The first list should contain a rectangle object for each of the sequences (this
        # will be the coloured box). We get this object from the bar objects that we are creating.
        # The second list is a list of the labels. This should be easier. We can just use the
        # most_abund_seq_names[:30] for these.
        # to grab the rectangles, I think its best to pick them up during the plotting and hold them in a list
        # outside of the subplot loop. We will need a holder for the objects to populate.

        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticks.html#matplotlib.axes.Axes.set_xticks
        axarr[axarr_index].set_yticks([1], minor=False)
        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html#matplotlib.axes.Axes.set_xticklabels
        axarr[axarr_index].set_yticklabels(['1'])
        axarr[axarr_index].set_xlabel(env_type)



        axarr_index += 1

    #https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots/26892326
    f.text(0, 0.5, 'ITS2 sequence relative abundance', va='center', rotation='vertical')
    ordered_rectange_list = [sub_list[0] for sub_list in legend_rectangle_holder]
    f.legend(ordered_rectange_list, top_n_seqs_to_get, loc='lower center')
    plt.tight_layout()
    f.savefig('diversity_bars.svg')
    f.show()

def create_diverstiy_figs_all_dates():
    ''' This is a development from the above code. Here in this code we will plot all of the date replicates
    for the env_samples as these could in theory aid a little bit in the diversity comparison.
    I will put these on the same axis as the samples of the same type with a gap or so inbetween.
    This will be the code to create the diversity figure for Gabi's paper.
        We have used symportal to do the QC work on the samples.
        As such we have tab delimited outputs that we can read into pandas dataframes
        We also have the info list from Gabi that contains what all of the samples are
        Esentially what we want to do is quite simple, however, we will need to be careful with the colouring
        As there are going to be a shed ton of sequences we will probably want to order the sequences found some how
        and colour those, then maybe colour sequences found at below 5% to grey. Something like that.
        I have a load of code from previous work when we were trying to look for actual SymPortal ITS2 type profiles
        within the env samples but we are no longer doing that. Still some of the code below might be useful.
        Actually I'm going to drop it into a junk code def and start from scratch here.
    '''

    try:
        sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
        QC_info_df= pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
        info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))
    except:
        # read in the SymPortal output
        sp_output_df = pd.read_csv('131_142.DIVs.absolute.txt', sep='\t', lineterminator='\n')

        # The SP output contains the QC info columns between the DIVs and the no_name ITS2 columns.
        # lets put the QC info columns into a seperate df.
        QC_info_df = sp_output_df[['Samples','raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique']]

        # now lets drop the QC columns from the SP output df and also drop the clade summation columns
        # we will be left with just clumns for each one of the sequences found in the samples
        sp_output_df.drop(columns=['noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                                'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                                'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique'
                                   ]
                          , inplace=True)

        # read in the info file
        info_df = pd.read_csv('info_300718.csv')

        # we no longer drop the dates
        # drop the rows that aren't from the 22.08.2016 data
        info_df = info_df[info_df['coral genus'] != 'Sponge']
        # need to use the ampersand rather than 'and': https://stackoverflow.com/questions/36921951/truth-value-of-a-series-is-ambiguous-use-a-empty-a-bool-a-item-a-any-o
        info_df = info_df[(info_df['Sample no.'] != 'negative extration') & (info_df['Sample no.'] != 'negative pcr') & (info_df['Sample no.'] != 'DIV_accession')]

        # Now we need to link the SP output to the sample names in the excel. Annoyingly they are formatted
        # slightly differently so we can't make a direct comparison.
        # easiest way to link them is to see if the first part of the SP name is the same as the first part
        # of the 'sequence file' in the meta info
        # when doing this we can also drop the SP info for those samples that won't be used i.e. those that
        # aren't now in the info_df

        # firstly rename the colums so that they are 'sample_name' in all of the dfs
        QC_info_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        sp_output_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        info_df.rename(index=str, columns={'Sample Name': 'sample_name'}, inplace=True)

        indices_to_drop = []
        for sp_index in sp_output_df.index.values.tolist():
            # keep track of whether the sp_index was found in the info table
            # if it wasn't then it should be dropped
            sys.stdout.write('\rsp_index: {}'.format(sp_index))
            found = False
            for info_index in info_df.index.values.tolist():
                if sp_output_df.loc[sp_index, 'sample_name'].split('_')[0] == info_df.loc[info_index, 'Sequence file'].split('_')[0]:
                    found = True
                    # then these are a related set of rows and we should make the sample_names the same
                    sp_output_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']
                    QC_info_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']


            if not found:
                indices_to_drop.append(sp_index)

        # drop the rows from the SP output tables that aren't going to be used
        sp_output_df.drop(inplace=True, index=indices_to_drop)
        QC_info_df.drop(inplace=True, index=indices_to_drop)

        # let's sort out the 'environ type' column in the info_df
        # currently it is a bit of a mess
        for index in info_df.index.values.tolist():
            if 'coral' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'coral'
            elif 'seawater' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sea_water'
            elif 'mucus' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'mucus'
            elif 'SA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_close'
            elif 'SB' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_far'
            elif 'TA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'turf'

        # now clean up the df indices
        info_df.index = range(len(info_df))
        sp_output_df.index = range(len(sp_output_df))
        QC_info_df.index = range(len(QC_info_df))

        # make the sample_name column the index for each of the datasets
        info_df.set_index('sample_name', inplace=True)
        sp_output_df.set_index('sample_name', inplace=True)
        QC_info_df.set_index('sample_name', inplace=True)

        # pickle the out put and put a check in place to see if we need to do the above
        pickle.dump(sp_output_df, open('sp_output_df_all_dates.pickle', 'wb'))
        pickle.dump(QC_info_df, open('QC_info_df_all_dates.pickle', 'wb'))
        pickle.dump(info_df, open('info_df_all_dates.pickle', 'wb'))



    # so we want to plot the ITS2 sequence diversity in each of the samples as bar charts
    # We are going to have a huge diversity of sequences to deal with so I think something along the lines
    # of plotting the top n most abundant sequences. The term 'most abundant' should be considered carefully here
    # I think it will be best if we work on a sample by sample basis. i.e. we pick the n sequences that have the
    # highest representation in any one sample. So for example what we are not doing is seeing how many times
    # C3 was sequenced across all of the samples, and finding that it is a lot and therefore plotting it.
    # we are looking in each of the samples and seeing the highest proportion it is found at in any one sample.
    # This way we should have the best chance of having a coloured representation for each sample's most abundant
    # sequence.

    # to start lets go sample by sample and see what the highest prop for each seq is.

    # dict to hold info on which sample and what the proportion is for each sequence
    # key = sequence name, value = tup ( sample name, relative abundance)
    try:
        seq_rel_abund_calculator_dict = pickle.load(open('seq_rel_abund_calculator_dict_all_dates.pickle', 'rb'))
    except:
        seq_rel_abund_calculator_dict = {}
        for sample_index in sp_output_df.index.values.tolist():
            sys.stdout.write('\nGeting rel seq abundances from {}\n'.format(sample_index))
            # temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum(axis='index'))
            temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum())
            for seq_name in temp_prop_array.keys():
                sys.stdout.write('\rseq: {}'.format(seq_name))
                val = temp_prop_array[seq_name]
                if val != 0:  # if the sequences was found in the sample
                    # If the sequence is already in the dict
                    if seq_name in seq_rel_abund_calculator_dict.keys():
                        # check to seee if the rel abundance is larger than the one already logged
                        if val > seq_rel_abund_calculator_dict[seq_name][1]:
                            seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
                    # if we haven't logged for this sequence yet, then add this as the first log
                    else:
                        seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
        pickle.dump(seq_rel_abund_calculator_dict, open('seq_rel_abund_calculator_dict_all_dates.pickle', 'wb'))



    # here we have a dict that contains the largest rel_abundances per sample for each of the seqs
    # now we can sort this to look at the top ? n sequences to start with (I'm not sure how the colouring will
    # look like so lets just start with 30 and see where we get to)
    sorted_list = sorted(seq_rel_abund_calculator_dict.items(), key = lambda x: x[1][1], reverse=True)
    most_abund_seq_names = [tup[0] for tup in sorted_list]

    # from the above sorted list we can then plot these sequences with colour and all others with grey scale
    # lets make a coulour dictionary for the most common types
    colour_list = get_colour_list()
    colour_dict = {}
    num_coloured_seqs = 30

    # we will also need a grey palette for those sequences that are not the ones being annotated
    grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    # make the dict
    for i in range(len(most_abund_seq_names)):
        if i < num_coloured_seqs:
            colour_dict[most_abund_seq_names[i]] = colour_list[i]
        else:
            colour_dict[most_abund_seq_names[i]] = grey_palette[i % 6]

    # add the 'low' key and assign it to grey for later on
    # the low category will be created later on to hold the grouped abundances of sequences in samples
    # below a certain rel abund cutoff
    colour_dict['low'] = '#D0CFD4'

    # for plotting we will also need to collect the 'rectangle' objects from the plotting process to use to make
    # the legend which we will display on the plot right at the end
    # this is going to be a little tricky to collect as not every sample type group that we are plotting
    # is going to have all of the top n sequences. So, we will have to pick up rectangles when they come up in
    # the plotting.
    # to keep track of which sequences we need we have:
    top_n_seqs_to_get = most_abund_seq_names[:num_coloured_seqs]
    # to keep track of which sequences we have already collected we have:
    top_n_seqs_done = []
    #finally to collect the rectangles and store them in order despite collecting them out of order we have:
    legend_rectangle_holder = [[] for i in range(num_coloured_seqs)]


    # set up the plotting environment
    # one plot for coral, mucus, seawater, sediment, turf_algae
    # f, axarr = plt.subplots(4, 5)
    # # counter to reference which set of axes we are plotting on
    # axarr_index = 0

    # rather than work with the subplot standard layout we will use subplot2grid for our set up and have some hard
    # coded axes
    # ax1 = (plt.subplot2grid((6, 5), (0, 0), colspan=3), 'coral', '22.08.2016')
    # ax2 = (plt.subplot2grid((6, 5), (0, 3), colspan=2), 'coral', '06.09.2016')
    # ax3 = (plt.subplot2grid((6, 5), (1, 0), colspan=3), 'mucus', '22.08.2016')


    ax_list = [
        # coral
        (plt.subplot2grid((6, 70), (0, 0), colspan=50), 'coral', '22.08.2016'),
        (plt.subplot2grid((6, 70), (0, 53), colspan=17), 'coral', '06.09.2016'),
        # mucus
        (plt.subplot2grid((6, 70), (1, 0), colspan=49), 'mucus', '22.08.2016'),
        (plt.subplot2grid((6, 70), (1, 53), colspan=17), 'mucus', '06.09.2016'),
        # sea water
        (plt.subplot2grid((6, 70), (2, 0 ),  colspan=5), 'sea_water', '22.08.2016'),
        (plt.subplot2grid((6, 70), (2, 13),  colspan=3), 'sea_water', '31.08.2016'),
        (plt.subplot2grid((6, 70), (2, 22), colspan=2), 'sea_water', '01.09.2016'),
        (plt.subplot2grid((6, 70), (2, 29), colspan=5), 'sea_water', '06.09.2016'),
        (plt.subplot2grid((6, 70), (2, 40), colspan=5), 'sea_water', '02.10.2016'),
        # sed_close
        (plt.subplot2grid((6, 70), (3, 0 ), colspan=5), 'sed_close', '22.08.2016'),
        (plt.subplot2grid((6, 70), (3, 13), colspan=3), 'sed_close', '31.08.2016'),
        (plt.subplot2grid((6, 70), (3, 22), colspan=2), 'sed_close', '01.09.2016'),
        (plt.subplot2grid((6, 70), (3, 29), colspan=5), 'sed_close', '06.09.2016'),
        (plt.subplot2grid((6, 70), (3, 40), colspan=5), 'sed_close', '02.10.2016'),
        # sed_far,
        (plt.subplot2grid((6, 70), (4, 0 ), colspan=5), 'sed_far', '22.08.2016'),
        (plt.subplot2grid((6, 70), (4, 13), colspan=3), 'sed_far', '31.08.2016'),
        (plt.subplot2grid((6, 70), (4, 22), colspan=2), 'sed_far', '01.09.2016'),
        (plt.subplot2grid((6, 70), (4, 29), colspan=5), 'sed_far', '06.09.2016'),
        (plt.subplot2grid((6, 70), (4, 40), colspan=5), 'sed_far', '02.10.2016'),
        # turf,
        (plt.subplot2grid((6, 70), (5, 0 ), colspan= 10), 'turf', '22.08.2016'),
        (plt.subplot2grid((6, 70), (5, 13), colspan=6 ), 'turf', '31.08.2016'),
        (plt.subplot2grid((6, 70), (5, 22), colspan=4 ), 'turf', '01.09.2016'),
        (plt.subplot2grid((6, 70), (5, 29), colspan=8 ), 'turf', '06.09.2016'),
        (plt.subplot2grid((6, 70), (5, 40), colspan=9 ), 'turf', '02.10.2016')
    ]



    # now we need to get the actual plotting information for each of the samples.
    # we can do this sample type by sample type
    # we will create a local dataframe that will be a sub set of the main sp output dataframe but will
    # only contain the samples of the given sample type. It will also eventually only contain the sequence
    # information for the sample type in question. Normally I would make a 2D list to hold the plotting
    # information but in this case having all of the information in a dataframe is working really well and
    # this is how I will do it in future.

    # go environment type by environment type
    # env_type_list = ['coral', 'mucus', 'sea_water', 'sed_close', 'sed_far', 'turf']
    # env_type_list = ['sea_water', 'sed_close', 'sed_far', 'turf']
    # for env_type in env_type_list:
    #     # then go date by date
    #     date_list = ['22.08.2016', '31.08.2016', '01.09.2016', '06.09.2016', '02.10.2016']
    #     for date in date_list:
    for ax_item in ax_list:
        try:
            env_sp_output_df_prop = pickle.load(open('env_sp_output_df_prop_{}_{}.pickle'.format(ax_item[1], ax_item[2]), 'rb'))
            sorted_list_of_env_specific_seqs = pickle.load(open('sorted_list_of_env_specific_seqs_{}_{}.pickle'.format(ax_item[1], ax_item[2]), 'rb'))
        except:
            sys.stdout.write('\nGenerating plotting info for {} samples {}\n'.format(ax_item[1], ax_item[2]))
            # currently we have something like 4000 sequences to plot which is too may
            # I think it will be much easier if we group the sequences that are found below a certain
            # threshold. I think the best way to do this is to create a slice of the main df that
            # contain the information for the samples of the env_type only

            # get subset of the main dfs that contain only the env_type samples from the given date
            env_info_df = info_df[(info_df['environ type'] == ax_item[1]) & (info_df['date collected'] == ax_item[2])]

            # we need to drop any samples that have only zeros in their columns
            # https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame
            env_info_df = env_info_df[(env_info_df.T != 0).any()]


            env_sp_output_df = sp_output_df.loc[env_info_df.index.values.tolist()]

            # we need to drop any samples that have only zeros in their columns
            # https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame
            env_sp_output_df = env_sp_output_df[(env_sp_output_df.T != 0).any()]

            # append a 'low' columns to the env_sp_ouput_df and populate with 0s
            env_sp_output_df['low'] = 0
            # now make a proportions version of the df, rather than absolute counts
            env_sp_output_df_prop = env_sp_output_df[:].div(env_sp_output_df[:].sum(axis=1), axis=0)

            # now as we work our way through we can sum up the low sequences into this column
            # we can then check for 0 columns and drop these.

            # get a list of the sequences found in the collection of samples of the given type and order
            # them according to summed rel_abundance acorss all samples. This should be the order in which
            # the samples are plotted
            # at the same time we can get the info we need for plotting

            summed_seq_rel_abund_across_smpls_dict = {seq: 0 for seq in list(sp_output_df)}

            # we are also going to need to put the samples in an order that makes sense.
            # Ideally I would like to take the time to do a proper travelling salesman analysis
            # and I was thinking of putting in an ant colony system to do this but...
            # I'm not sure we've got time for that. So let's do something a little more simple which
            # should still be effective and look groovy
            # Lets sort according to the majority sequence. I.e. lets put the samples into groups that
            # are defined by what their majorty sequence is, then lets plot in order of those groups.
            # within the groups we will plot in the order of the most abundant rel abund first.
            # to get this information we will have a dict to collect it. The key of the dict will
            # be the majority sequence with the value being a list that contains tuples. One tuple
            # for each sample in the list which will contain the sample_name and the re abund of the maj seq

            # Dict for holding the sample_sorting info
            sample_sorting_info_dict = defaultdict(list)

            for sp_index in env_sp_output_df_prop.index.values.tolist():
                sys.stdout.write('\r{}'.format(sp_index))
                # we need to get the name of the most abundant sequence and its rel abund for each sample
                # for the sample_sorting_info_dict
                most_abund_seq_name = env_sp_output_df_prop.loc[sp_index].idxmax(axis='index')

                rel_abund_of_most_abund_seq = env_sp_output_df_prop.loc[sp_index, most_abund_seq_name]

                apples = 'asdf'
                sample_sorting_info_dict[most_abund_seq_name].append((sp_index, rel_abund_of_most_abund_seq))
                # Put its seq rel abundances ino the summed_seq_rel_abund_across... dict

                for non_zero_seq in env_sp_output_df_prop.loc[sp_index][env_sp_output_df_prop.loc[sp_index] > 0].index:  # not including the final 'low' columns (this will be zero)
                    val = env_sp_output_df_prop.loc[sp_index, non_zero_seq]
                    # be sure to count the value of this cell and using it in judging which are the most
                    # abundant sequences before we check whether to relegate it to the 'low' column
                    summed_seq_rel_abund_across_smpls_dict[non_zero_seq] += val
                    if val < 0.005:
                        env_sp_output_df_prop.loc[sp_index, 'low'] += val
                        env_sp_output_df_prop.loc[sp_index, non_zero_seq] = 0


            # here we can get a sorted sample list using the sample_sorting_info_dict
            sorted_sample_list = []
            # we want to work through the sample_sorting_info_dict by the longest lists first
            sorted_keys_for_sample_sort = \
                [tup[0] for tup in sorted(sample_sorting_info_dict.items(), key=lambda x: len(x[1]), reverse=True)]
            # for each of teh maj seq groups
            for sorted_key in sorted_keys_for_sample_sort:
                # now within each of these lists we want to order according to the rel_abundance of the sequences
                sorted_list_of_samples_of_group = [tup[0] for tup in sorted(sample_sorting_info_dict[sorted_key], key=lambda x: x[1], reverse=True)]
                sorted_sample_list.extend(sorted_list_of_samples_of_group)

            # now we should re-order the df so that it is in the sample order of sorted_sample_list
            env_sp_output_df_prop = env_sp_output_df_prop.reindex(sorted_sample_list)

            # here we have a dict that contains the abundances of the seqs for the coral samples
            # we will plot the coral samples' seqs in the order of the sequences in this dict
            sorted_list_of_env_specific_seqs_tup \
                = sorted(summed_seq_rel_abund_across_smpls_dict.items(), key=lambda x: x[1], reverse=True)
            sorted_list_of_env_specific_seqs = [tup[0] for tup in sorted_list_of_env_specific_seqs_tup]

            # Now we check for zero only columns and drop them
            # we also need to remove these columns from the sorted
            non_zero_cols = list(env_sp_output_df_prop.loc[:, (env_sp_output_df_prop != 0).any(axis=0)])
            sorted_list_of_env_specific_seqs = [seq for seq in sorted_list_of_env_specific_seqs if seq in non_zero_cols]

            # add the 'low' column only if the 'low' column is a non-zero column, i.e. used
            if 'low' in non_zero_cols:
                sorted_list_of_env_specific_seqs.append('low')

            # now drop the cols
            env_sp_output_df_prop = env_sp_output_df_prop[non_zero_cols]

            # we also have the plotting list which contains the info we will be plotting.

            # the plotting_list is currently organised in a different order to that of the sorted_list_of_env...
            # we need to change this order
            env_sp_output_df_prop = env_sp_output_df_prop[sorted_list_of_env_specific_seqs]

            # we now need to transpose this
            env_sp_output_df_prop = env_sp_output_df_prop.transpose()

            # pickle the info we've been plotting
            pickle.dump(env_sp_output_df_prop, open('env_sp_output_df_prop_{}_{}.pickle'.format(ax_item[1], ax_item[2]), 'wb'))
            pickle.dump(sorted_list_of_env_specific_seqs, open('sorted_list_of_env_specific_seqs_{}_{}.pickle'.format(ax_item[1], ax_item[2]), 'wb'))

        # here we finally have the plotting lists in the order of the sorted_list_of_env
        bottom = [0 for smp in list(env_sp_output_df_prop)]
        bar_list = []
        # for each sample
        plotting_indices = range(len(list(env_sp_output_df_prop)))

        # for each sequence
        list_of_seqs = env_sp_output_df_prop.index.values.tolist()
        for i in range(len(list_of_seqs)):

            sys.stdout.write('\r{}/{}'.format(i, len(list(env_sp_output_df_prop.iloc[:, 0]))))
            bar_list.append(ax_item[0].bar(plotting_indices, list(env_sp_output_df_prop.iloc[i]), 1, bottom,
                                                   color=colour_dict[sorted_list_of_env_specific_seqs[i]]))
            bottom = [L + M for L, M in zip(bottom, list(env_sp_output_df_prop.iloc[i]))]

            # check to see if the seq we are plotting is still in the top_n_seqs list
            # if it is in this list then we still need to grab a rectangle for plotting
            # the legend. Once we have grabbed the rectangle then we should remove the seq from the top_n_seqs list
            seq_name = list_of_seqs[i]

            if seq_name in top_n_seqs_to_get and seq_name not in top_n_seqs_done:
                # then this is a seq that we still need to get a rectangle for the legend
                legend_rectangle_holder[top_n_seqs_to_get.index(seq_name)].append(bar_list[-1][0])
                top_n_seqs_done.append(seq_name)

        # https://stackoverflow.com/questions/12998430/remove-xticks-in-a-matplotlib-plot
        ax_item[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        # https://stackoverflow.com/questions/15858192/how-to-set-xlim-and-ylim-for-a-subplot-in-matplotlib
        ax_item[0].set_xlim((-0.5, (len(list(env_sp_output_df_prop)) - 0.5)))
        # https://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib
        ax_item[0].spines['right'].set_visible(False)
        ax_item[0].spines['top'].set_visible(False)

        # To make the legend for this mo fo is going to be a little tricky. The legend argument basically takes
        # two lists. The first list should contain a rectangle object for each of the sequences (this
        # will be the coloured box). We get this object from the bar objects that we are creating.
        # The second list is a list of the labels. This should be easier. We can just use the
        # most_abund_seq_names[:30] for these.
        # to grab the rectangles, I think its best to pick them up during the plotting and hold them in a list
        # outside of the subplot loop. We will need a holder for the objects to populate.

        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticks.html#matplotlib.axes.Axes.set_xticks
        ax_item[0].set_yticks([1], minor=False)
        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html#matplotlib.axes.Axes.set_xticklabels
        ax_item[0].set_yticklabels(['1'])
        if ax_item[2] == '22.08.2016':
            ax_item[0].set_ylabel(ax_item[1])

        # only if we're on the final plot then we should add the date ticks
        if ax_item[1] == 'turf':
            ax_item[0].set_xlabel(ax_item[2])
        # also if we are on the mucus we should put the labels for the dates on the axis
        if ax_item[1] == 'mucus':
            ax_item[0].set_xlabel(ax_item[2])




    #https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots/26892326
    plt.text(0, 0.5, 'ITS2 sequence relative abundance', va='center', rotation='vertical')
    ordered_rectange_list = [sub_list[0] for sub_list in legend_rectangle_holder]
    plt.legend(ordered_rectange_list, top_n_seqs_to_get, loc='lower center')
    #plt.tight_layout()
    plt.savefig('env_smp_diversity_bars_all_dates.svg')
    plt.show()

def create_diverstiy_figs_all_dates_grouped_dates():
    '''
    This is yet another development in which we will group all of the dates together rather than sorting them
    individually.
    This is a development from the above code. Here in this code we will plot all of the date replicates
    for the env_samples as these could in theory aid a little bit in the diversity comparison.
    I will put these on the same axis as the samples of the same type with a gap or so inbetween.
    This will be the code to create the diversity figure for Gabi's paper.
        We have used symportal to do the QC work on the samples.
        As such we have tab delimited outputs that we can read into pandas dataframes
        We also have the info list from Gabi that contains what all of the samples are
        Esentially what we want to do is quite simple, however, we will need to be careful with the colouring
        As there are going to be a shed ton of sequences we will probably want to order the sequences found some how
        and colour those, then maybe colour sequences found at below 5% to grey. Something like that.
        I have a load of code from previous work when we were trying to look for actual SymPortal ITS2 type profiles
        within the env samples but we are no longer doing that. Still some of the code below might be useful.
        Actually I'm going to drop it into a junk code def and start from scratch here.
    '''

    try:
        sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
        QC_info_df= pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
        info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))
    except:
        # read in the SymPortal output
        sp_output_df = pd.read_csv('131_142.DIVs.absolute.txt', sep='\t', lineterminator='\n')

        # The SP output contains the QC info columns between the DIVs and the no_name ITS2 columns.
        # lets put the QC info columns into a seperate df.
        QC_info_df = sp_output_df[['Samples','raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique']]

        # now lets drop the QC columns from the SP output df and also drop the clade summation columns
        # we will be left with just clumns for each one of the sequences found in the samples
        sp_output_df.drop(columns=['noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                                'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                                'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique'
                                   ]
                          , inplace=True)

        # read in the info file
        info_df = pd.read_csv('info_300718.csv')

        # we no longer drop the dates
        # drop the rows that aren't from the 22.08.2016 data
        info_df = info_df[info_df['coral genus'] != 'Sponge']
        # need to use the ampersand rather than 'and': https://stackoverflow.com/questions/36921951/truth-value-of-a-series-is-ambiguous-use-a-empty-a-bool-a-item-a-any-o
        info_df = info_df[(info_df['Sample no.'] != 'negative extration') & (info_df['Sample no.'] != 'negative pcr') & (info_df['Sample no.'] != 'DIV_accession')]

        # Now we need to link the SP output to the sample names in the excel. Annoyingly they are formatted
        # slightly differently so we can't make a direct comparison.
        # easiest way to link them is to see if the first part of the SP name is the same as the first part
        # of the 'sequence file' in the meta info
        # when doing this we can also drop the SP info for those samples that won't be used i.e. those that
        # aren't now in the info_df

        # firstly rename the colums so that they are 'sample_name' in all of the dfs
        QC_info_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        sp_output_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        info_df.rename(index=str, columns={'Sample Name': 'sample_name'}, inplace=True)

        indices_to_drop = []
        for sp_index in sp_output_df.index.values.tolist():
            # keep track of whether the sp_index was found in the info table
            # if it wasn't then it should be dropped
            sys.stdout.write('\rsp_index: {}'.format(sp_index))
            found = False
            for info_index in info_df.index.values.tolist():
                if sp_output_df.loc[sp_index, 'sample_name'].split('_')[0] == info_df.loc[info_index, 'Sequence file'].split('_')[0]:
                    found = True
                    # then these are a related set of rows and we should make the sample_names the same
                    sp_output_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']
                    QC_info_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']


            if not found:
                indices_to_drop.append(sp_index)

        # drop the rows from the SP output tables that aren't going to be used
        sp_output_df.drop(inplace=True, index=indices_to_drop)
        QC_info_df.drop(inplace=True, index=indices_to_drop)

        # let's sort out the 'environ type' column in the info_df
        # currently it is a bit of a mess
        for index in info_df.index.values.tolist():
            if 'coral' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'coral'
            elif 'seawater' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sea_water'
            elif 'mucus' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'mucus'
            elif 'SA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_close'
            elif 'SB' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_far'
            elif 'TA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'turf'

        # now clean up the df indices
        info_df.index = range(len(info_df))
        sp_output_df.index = range(len(sp_output_df))
        QC_info_df.index = range(len(QC_info_df))

        # make the sample_name column the index for each of the datasets
        info_df.set_index('sample_name', inplace=True)
        sp_output_df.set_index('sample_name', inplace=True)
        QC_info_df.set_index('sample_name', inplace=True)

        # pickle the out put and put a check in place to see if we need to do the above
        pickle.dump(sp_output_df, open('sp_output_df_all_dates.pickle', 'wb'))
        pickle.dump(QC_info_df, open('QC_info_df_all_dates.pickle', 'wb'))
        pickle.dump(info_df, open('info_df_all_dates.pickle', 'wb'))



    # so we want to plot the ITS2 sequence diversity in each of the samples as bar charts
    # We are going to have a huge diversity of sequences to deal with so I think something along the lines
    # of plotting the top n most abundant sequences. The term 'most abundant' should be considered carefully here
    # I think it will be best if we work on a sample by sample basis. i.e. we pick the n sequences that have the
    # highest representation in any one sample. So for example what we are not doing is seeing how many times
    # C3 was sequenced across all of the samples, and finding that it is a lot and therefore plotting it.
    # we are looking in each of the samples and seeing the highest proportion it is found at in any one sample.
    # This way we should have the best chance of having a coloured representation for each sample's most abundant
    # sequence.

    # to start lets go sample by sample and see what the highest prop for each seq is.

    # dict to hold info on which sample and what the proportion is for each sequence
    # key = sequence name, value = tup ( sample name, relative abundance)
    try:
        seq_rel_abund_calculator_dict = pickle.load(open('seq_rel_abund_calculator_dict_all_dates.pickle', 'rb'))
    except:
        seq_rel_abund_calculator_dict = {}
        for sample_index in sp_output_df.index.values.tolist():
            sys.stdout.write('\nGeting rel seq abundances from {}\n'.format(sample_index))
            # temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum(axis='index'))
            temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum())
            for seq_name in temp_prop_array.keys():
                sys.stdout.write('\rseq: {}'.format(seq_name))
                val = temp_prop_array[seq_name]
                if val != 0:  # if the sequences was found in the sample
                    # If the sequence is already in the dict
                    if seq_name in seq_rel_abund_calculator_dict.keys():
                        # check to seee if the rel abundance is larger than the one already logged
                        if val > seq_rel_abund_calculator_dict[seq_name][1]:
                            seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
                    # if we haven't logged for this sequence yet, then add this as the first log
                    else:
                        seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
        pickle.dump(seq_rel_abund_calculator_dict, open('seq_rel_abund_calculator_dict_all_dates.pickle', 'wb'))



    # here we have a dict that contains the largest rel_abundances per sample for each of the seqs
    # now we can sort this to look at the top ? n sequences to start with (I'm not sure how the colouring will
    # look like so lets just start with 30 and see where we get to)
    sorted_list = sorted(seq_rel_abund_calculator_dict.items(), key = lambda x: x[1][1], reverse=True)
    most_abund_seq_names = [tup[0] for tup in sorted_list]

    # from the above sorted list we can then plot these sequences with colour and all others with grey scale
    # lets make a coulour dictionary for the most common types
    colour_list = get_colour_list()
    colour_dict = {}
    num_coloured_seqs = 30

    # we will also need a grey palette for those sequences that are not the ones being annotated
    grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    # make the dict
    for i in range(len(most_abund_seq_names)):
        if i < num_coloured_seqs:
            colour_dict[most_abund_seq_names[i]] = colour_list[i]
        else:
            colour_dict[most_abund_seq_names[i]] = grey_palette[i % 6]

    # add the 'low' key and assign it to grey for later on
    # the low category will be created later on to hold the grouped abundances of sequences in samples
    # below a certain rel abund cutoff
    colour_dict['low'] = '#D0CFD4'

    # for plotting we will also need to collect the 'rectangle' objects from the plotting process to use to make
    # the legend which we will display on the plot right at the end
    # this is going to be a little tricky to collect as not every sample type group that we are plotting
    # is going to have all of the top n sequences. So, we will have to pick up rectangles when they come up in
    # the plotting.
    # to keep track of which sequences we need we have:
    top_n_seqs_to_get = most_abund_seq_names[:num_coloured_seqs]
    # to keep track of which sequences we have already collected we have:
    top_n_seqs_done = []
    #finally to collect the rectangles and store them in order despite collecting them out of order we have:
    legend_rectangle_holder = [[] for i in range(num_coloured_seqs)]


    # set up the plotting environment


    # rather than work with the subplot standard layout we will use subplot2grid for our set up and have some hard
    # coded axes

    this_fig = plt.figure(figsize=(10, 8))

    ax_list = [
        # coral
        (plt.subplot2grid((6, 67), (0, 0), colspan=67), 'coral'),

        # mucus
        (plt.subplot2grid((6, 67), (1, 0), colspan=66), 'mucus'),

        # sea water
        (plt.subplot2grid((6, 67), (2, 0 ),  colspan=20), 'sea_water'),

        # sed_close
        (plt.subplot2grid((6, 67), (3, 0 ), colspan=20), 'sed_close'),

        # sed_far,
        (plt.subplot2grid((6, 67), (4, 0 ), colspan=20), 'sed_far'),

        # turf,
        (plt.subplot2grid((6, 67), (5, 0 ), colspan= 37), 'turf'),

    ]





    # now we need to get the actual plotting information for each of the samples.
    # we can do this sample type by sample type
    # we will create a local dataframe that will be a sub set of the main sp output dataframe but will
    # only contain the samples of the given sample type. It will also eventually only contain the sequence
    # information for the sample type in question. Normally I would make a 2D list to hold the plotting
    # information but in this case having all of the information in a dataframe is working really well and
    # this is how I will do it in future.

    # go environment type by environment type
    # env_type_list = ['coral', 'mucus', 'sea_water', 'sed_close', 'sed_far', 'turf']
    # env_type_list = ['sea_water', 'sed_close', 'sed_far', 'turf']
    # for env_type in env_type_list:
    #     # then go date by date
    #     date_list = ['22.08.2016', '31.08.2016', '01.09.2016', '06.09.2016', '02.10.2016']
    #     for date in date_list:
    for ax_item in ax_list:
        try:
            env_sp_output_df_prop = pickle.load(open('env_sp_output_df_prop_{}_grouped.pickle'.format(ax_item[1]), 'rb'))
            sorted_list_of_env_specific_seqs = pickle.load(open('sorted_list_of_env_specific_seqs_{}_grouped.pickle'.format(ax_item[1]), 'rb'))
        except:
            sys.stdout.write('\nGenerating plotting info for samples {}\n'.format(ax_item[1]))
            # currently we have something like 4000 sequences to plot which is too may
            # I think it will be much easier if we group the sequences that are found below a certain
            # threshold. I think the best way to do this is to create a slice of the main df that
            # contain the information for the samples of the env_type only

            # get subset of the main dfs that contain only the env_type samples from the given date
            env_info_df = info_df[(info_df['environ type'] == ax_item[1])]

            # we need to drop any samples that have only zeros in their columns
            # https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame
            env_info_df = env_info_df[(env_info_df.T != 0).any()]


            env_sp_output_df = sp_output_df.loc[env_info_df.index.values.tolist()]

            # we need to drop any samples that have only zeros in their columns
            # https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame
            env_sp_output_df = env_sp_output_df[(env_sp_output_df.T != 0).any()]

            # append a 'low' columns to the env_sp_ouput_df and populate with 0s
            env_sp_output_df['low'] = 0
            # now make a proportions version of the df, rather than absolute counts
            env_sp_output_df_prop = env_sp_output_df[:].div(env_sp_output_df[:].sum(axis=1), axis=0)

            # now as we work our way through we can sum up the low sequences into this column
            # we can then check for 0 columns and drop these.

            # get a list of the sequences found in the collection of samples of the given type and order
            # them according to summed rel_abundance acorss all samples. This should be the order in which
            # the samples are plotted
            # at the same time we can get the info we need for plotting

            summed_seq_rel_abund_across_smpls_dict = {seq: 0 for seq in list(sp_output_df)}

            # we are also going to need to put the samples in an order that makes sense.
            # Ideally I would like to take the time to do a proper travelling salesman analysis
            # and I was thinking of putting in an ant colony system to do this but...
            # I'm not sure we've got time for that. So let's do something a little more simple which
            # should still be effective and look groovy
            # Lets sort according to the majority sequence. I.e. lets put the samples into groups that
            # are defined by what their majorty sequence is, then lets plot in order of those groups.
            # within the groups we will plot in the order of the most abundant rel abund first.
            # to get this information we will have a dict to collect it. The key of the dict will
            # be the majority sequence with the value being a list that contains tuples. One tuple
            # for each sample in the list which will contain the sample_name and the re abund of the maj seq

            # Dict for holding the sample_sorting info
            sample_sorting_info_dict = defaultdict(list)

            for sp_index in env_sp_output_df_prop.index.values.tolist():
                sys.stdout.write('\r{}'.format(sp_index))
                # we need to get the name of the most abundant sequence and its rel abund for each sample
                # for the sample_sorting_info_dict
                most_abund_seq_name = env_sp_output_df_prop.loc[sp_index].idxmax(axis='index')

                rel_abund_of_most_abund_seq = env_sp_output_df_prop.loc[sp_index, most_abund_seq_name]

                apples = 'asdf'
                sample_sorting_info_dict[most_abund_seq_name].append((sp_index, rel_abund_of_most_abund_seq))
                # Put its seq rel abundances ino the summed_seq_rel_abund_across... dict

                for non_zero_seq in env_sp_output_df_prop.loc[sp_index][env_sp_output_df_prop.loc[sp_index] > 0].index:  # not including the final 'low' columns (this will be zero)
                    val = env_sp_output_df_prop.loc[sp_index, non_zero_seq]
                    # be sure to count the value of this cell and using it in judging which are the most
                    # abundant sequences before we check whether to relegate it to the 'low' column
                    summed_seq_rel_abund_across_smpls_dict[non_zero_seq] += val
                    if val < 0.005:
                        env_sp_output_df_prop.loc[sp_index, 'low'] += val
                        env_sp_output_df_prop.loc[sp_index, non_zero_seq] = 0


            # here we can get a sorted sample list using the sample_sorting_info_dict
            sorted_sample_list = []
            # we want to work through the sample_sorting_info_dict by the longest lists first
            sorted_keys_for_sample_sort = \
                [tup[0] for tup in sorted(sample_sorting_info_dict.items(), key=lambda x: len(x[1]), reverse=True)]
            # for each of teh maj seq groups
            for sorted_key in sorted_keys_for_sample_sort:
                # now within each of these lists we want to order according to the rel_abundance of the sequences
                sorted_list_of_samples_of_group = [tup[0] for tup in sorted(sample_sorting_info_dict[sorted_key], key=lambda x: x[1], reverse=True)]
                sorted_sample_list.extend(sorted_list_of_samples_of_group)

            # now we should re-order the df so that it is in the sample order of sorted_sample_list
            env_sp_output_df_prop = env_sp_output_df_prop.reindex(sorted_sample_list)

            # here we have a dict that contains the abundances of the seqs for the coral samples
            # we will plot the coral samples' seqs in the order of the sequences in this dict
            sorted_list_of_env_specific_seqs_tup \
                = sorted(summed_seq_rel_abund_across_smpls_dict.items(), key=lambda x: x[1], reverse=True)
            sorted_list_of_env_specific_seqs = [tup[0] for tup in sorted_list_of_env_specific_seqs_tup]

            # Now we check for zero only columns and drop them
            # we also need to remove these columns from the sorted
            non_zero_cols = list(env_sp_output_df_prop.loc[:, (env_sp_output_df_prop != 0).any(axis=0)])
            sorted_list_of_env_specific_seqs = [seq for seq in sorted_list_of_env_specific_seqs if seq in non_zero_cols]

            # add the 'low' column only if the 'low' column is a non-zero column, i.e. used
            if 'low' in non_zero_cols:
                sorted_list_of_env_specific_seqs.append('low')

            # now drop the cols
            env_sp_output_df_prop = env_sp_output_df_prop[non_zero_cols]

            # we also have the plotting list which contains the info we will be plotting.

            # the plotting_list is currently organised in a different order to that of the sorted_list_of_env...
            # we need to change this order
            env_sp_output_df_prop = env_sp_output_df_prop[sorted_list_of_env_specific_seqs]

            # we now need to transpose this
            env_sp_output_df_prop = env_sp_output_df_prop.transpose()

            # pickle the info we've been plotting
            pickle.dump(env_sp_output_df_prop, open('env_sp_output_df_prop_{}_grouped.pickle'.format(ax_item[1]), 'wb'))
            pickle.dump(sorted_list_of_env_specific_seqs, open('sorted_list_of_env_specific_seqs_{}_grouped.pickle'.format(ax_item[1]), 'wb'))

        # here we finally have the plotting lists in the order of the sorted_list_of_env
        bottom = [0 for smp in list(env_sp_output_df_prop)]
        bar_list = []
        # for each sample
        plotting_indices = range(len(list(env_sp_output_df_prop)))

        # for each sequence
        list_of_seqs = env_sp_output_df_prop.index.values.tolist()
        for i in range(len(list_of_seqs)):

            sys.stdout.write('\r{}/{}'.format(i, len(list(env_sp_output_df_prop.iloc[:, 0]))))
            bar_list.append(ax_item[0].bar(plotting_indices, list(env_sp_output_df_prop.iloc[i]), 1, bottom,
                                                   color=colour_dict[sorted_list_of_env_specific_seqs[i]]))
            bottom = [L + M for L, M in zip(bottom, list(env_sp_output_df_prop.iloc[i]))]

            # check to see if the seq we are plotting is still in the top_n_seqs list
            # if it is in this list then we still need to grab a rectangle for plotting
            # the legend. Once we have grabbed the rectangle then we should remove the seq from the top_n_seqs list
            seq_name = list_of_seqs[i]

            if seq_name in top_n_seqs_to_get and seq_name not in top_n_seqs_done:
                # then this is a seq that we still need to get a rectangle for the legend
                legend_rectangle_holder[top_n_seqs_to_get.index(seq_name)].append(bar_list[-1][0])
                top_n_seqs_done.append(seq_name)

        # https://stackoverflow.com/questions/12998430/remove-xticks-in-a-matplotlib-plot
        ax_item[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        # https://stackoverflow.com/questions/15858192/how-to-set-xlim-and-ylim-for-a-subplot-in-matplotlib
        ax_item[0].set_xlim((-0.5, (len(list(env_sp_output_df_prop)) - 0.5)))
        # https://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib
        ax_item[0].spines['right'].set_visible(False)
        ax_item[0].spines['top'].set_visible(False)

        # To make the legend for this mo fo is going to be a little tricky. The legend argument basically takes
        # two lists. The first list should contain a rectangle object for each of the sequences (this
        # will be the coloured box). We get this object from the bar objects that we are creating.
        # The second list is a list of the labels. This should be easier. We can just use the
        # most_abund_seq_names[:30] for these.
        # to grab the rectangles, I think its best to pick them up during the plotting and hold them in a list
        # outside of the subplot loop. We will need a holder for the objects to populate.

        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticks.html#matplotlib.axes.Axes.set_xticks
        ax_item[0].set_yticks([1], minor=False)
        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html#matplotlib.axes.Axes.set_xticklabels
        ax_item[0].set_yticklabels(['1'])
        # if ax_item[2] == '22.08.2016':
        ax_item[0].set_ylabel(ax_item[1])

        # # only if we're on the final plot then we should add the date ticks
        # if ax_item[1] == 'turf':
        #     ax_item[0].set_xlabel(ax_item[2])
        # # also if we are on the mucus we should put the labels for the dates on the axis
        # if ax_item[1] == 'mucus':
        #     ax_item[0].set_xlabel(ax_item[2])




    #https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots/26892326
    plt.text(0, 0.5, 'ITS2 sequence relative abundance', va='center', rotation='vertical')
    ordered_rectange_list = [sub_list[0] for sub_list in legend_rectangle_holder]
    plt.legend(ordered_rectange_list, top_n_seqs_to_get, loc='lower center')
    #plt.tight_layout()
    plt.savefig('env_smp_diversity_bars_all_dates_grouped.svg')
    plt.savefig('env_smp_diversity_bars_all_dates_grouped.png')
    # plt.show()

def create_diverstiy_figs_all_dates_grouped_dates_sed_grouped():
    '''
    This is yet another development in which we will group all of the dates together rather than sorting them
    individually.
    This is a development from the above code. Here in this code we will plot all of the date replicates
    for the env_samples as these could in theory aid a little bit in the diversity comparison.
    I will put these on the same axis as the samples of the same type with a gap or so inbetween.
    This will be the code to create the diversity figure for Gabi's paper.
        We have used symportal to do the QC work on the samples.
        As such we have tab delimited outputs that we can read into pandas dataframes
        We also have the info list from Gabi that contains what all of the samples are
        Esentially what we want to do is quite simple, however, we will need to be careful with the colouring
        As there are going to be a shed ton of sequences we will probably want to order the sequences found some how
        and colour those, then maybe colour sequences found at below 5% to grey. Something like that.
        I have a load of code from previous work when we were trying to look for actual SymPortal ITS2 type profiles
        within the env samples but we are no longer doing that. Still some of the code below might be useful.
        Actually I'm going to drop it into a junk code def and start from scratch here.
    '''

    try:
        sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
        QC_info_df= pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
        info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))
    except:
        # read in the SymPortal output
        sp_output_df = pd.read_csv('131_142.DIVs.absolute.txt', sep='\t', lineterminator='\n')

        # The SP output contains the QC info columns between the DIVs and the no_name ITS2 columns.
        # lets put the QC info columns into a seperate df.
        QC_info_df = sp_output_df[['Samples','raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique']]

        # now lets drop the QC columns from the SP output df and also drop the clade summation columns
        # we will be left with just clumns for each one of the sequences found in the samples
        sp_output_df.drop(columns=['noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                                'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                                'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                                'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                                'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
                                   'size_screening_violation_absolute', 'size_screening_violation_unique',
                                   'post_med_absolute', 'post_med_unique'
                                   ]
                          , inplace=True)

        # read in the info file
        info_df = pd.read_csv('info_300718.csv')

        # we no longer drop the dates
        # drop the rows that aren't from the 22.08.2016 data
        info_df = info_df[info_df['coral genus'] != 'Sponge']
        # need to use the ampersand rather than 'and': https://stackoverflow.com/questions/36921951/truth-value-of-a-series-is-ambiguous-use-a-empty-a-bool-a-item-a-any-o
        info_df = info_df[(info_df['Sample no.'] != 'negative extration') & (info_df['Sample no.'] != 'negative pcr') & (info_df['Sample no.'] != 'DIV_accession')]

        # Now we need to link the SP output to the sample names in the excel. Annoyingly they are formatted
        # slightly differently so we can't make a direct comparison.
        # easiest way to link them is to see if the first part of the SP name is the same as the first part
        # of the 'sequence file' in the meta info
        # when doing this we can also drop the SP info for those samples that won't be used i.e. those that
        # aren't now in the info_df

        # firstly rename the colums so that they are 'sample_name' in all of the dfs
        QC_info_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        sp_output_df.rename(index=str, columns={'Samples': 'sample_name'}, inplace=True)
        info_df.rename(index=str, columns={'Sample Name': 'sample_name'}, inplace=True)

        indices_to_drop = []
        for sp_index in sp_output_df.index.values.tolist():
            # keep track of whether the sp_index was found in the info table
            # if it wasn't then it should be dropped
            sys.stdout.write('\rsp_index: {}'.format(sp_index))
            found = False
            for info_index in info_df.index.values.tolist():
                if sp_output_df.loc[sp_index, 'sample_name'].split('_')[0] == info_df.loc[info_index, 'Sequence file'].split('_')[0]:
                    found = True
                    # then these are a related set of rows and we should make the sample_names the same
                    sp_output_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']
                    QC_info_df.loc[sp_index, 'sample_name'] = info_df.loc[info_index, 'sample_name']


            if not found:
                indices_to_drop.append(sp_index)

        # drop the rows from the SP output tables that aren't going to be used
        sp_output_df.drop(inplace=True, index=indices_to_drop)
        QC_info_df.drop(inplace=True, index=indices_to_drop)

        # let's sort out the 'environ type' column in the info_df
        # currently it is a bit of a mess
        for index in info_df.index.values.tolist():
            if 'coral' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'coral'
            elif 'seawater' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sea_water'
            elif 'mucus' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'mucus'
            elif 'SA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_close'
            elif 'SB' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'sed_far'
            elif 'TA' in info_df.loc[index, 'Sample no.']:
                info_df.loc[index, 'environ type'] = 'turf'

        # now clean up the df indices
        info_df.index = range(len(info_df))
        sp_output_df.index = range(len(sp_output_df))
        QC_info_df.index = range(len(QC_info_df))

        # make the sample_name column the index for each of the datasets
        info_df.set_index('sample_name', inplace=True)
        sp_output_df.set_index('sample_name', inplace=True)
        QC_info_df.set_index('sample_name', inplace=True)

        # pickle the out put and put a check in place to see if we need to do the above
        pickle.dump(sp_output_df, open('sp_output_df_all_dates.pickle', 'wb'))
        pickle.dump(QC_info_df, open('QC_info_df_all_dates.pickle', 'wb'))
        pickle.dump(info_df, open('info_df_all_dates.pickle', 'wb'))


    info_df = info_df.applymap(lambda x: 'sed' if "sed" in str(x) else x)

    # so we want to plot the ITS2 sequence diversity in each of the samples as bar charts
    # We are going to have a huge diversity of sequences to deal with so I think something along the lines
    # of plotting the top n most abundant sequences. The term 'most abundant' should be considered carefully here
    # I think it will be best if we work on a sample by sample basis. i.e. we pick the n sequences that have the
    # highest representation in any one sample. So for example what we are not doing is seeing how many times
    # C3 was sequenced across all of the samples, and finding that it is a lot and therefore plotting it.
    # we are looking in each of the samples and seeing the highest proportion it is found at in any one sample.
    # This way we should have the best chance of having a coloured representation for each sample's most abundant
    # sequence.

    # to start lets go sample by sample and see what the highest prop for each seq is.

    # dict to hold info on which sample and what the proportion is for each sequence
    # key = sequence name, value = tup ( sample name, relative abundance)
    try:
        seq_rel_abund_calculator_dict = pickle.load(open('seq_rel_abund_calculator_dict_all_dates.pickle', 'rb'))
    except:
        seq_rel_abund_calculator_dict = {}
        for sample_index in sp_output_df.index.values.tolist():
            sys.stdout.write('\nGeting rel seq abundances from {}\n'.format(sample_index))
            # temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum(axis='index'))
            temp_prop_array = sp_output_df.loc[sample_index].div(sp_output_df.loc[sample_index].sum())
            for seq_name in temp_prop_array.keys():
                sys.stdout.write('\rseq: {}'.format(seq_name))
                val = temp_prop_array[seq_name]
                if val != 0:  # if the sequences was found in the sample
                    # If the sequence is already in the dict
                    if seq_name in seq_rel_abund_calculator_dict.keys():
                        # check to seee if the rel abundance is larger than the one already logged
                        if val > seq_rel_abund_calculator_dict[seq_name][1]:
                            seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
                    # if we haven't logged for this sequence yet, then add this as the first log
                    else:
                        seq_rel_abund_calculator_dict[seq_name] = (sample_index, val)
        pickle.dump(seq_rel_abund_calculator_dict, open('seq_rel_abund_calculator_dict_all_dates.pickle', 'wb'))



    # here we have a dict that contains the largest rel_abundances per sample for each of the seqs
    # now we can sort this to look at the top ? n sequences to start with (I'm not sure how the colouring will
    # look like so lets just start with 30 and see where we get to)
    sorted_list = sorted(seq_rel_abund_calculator_dict.items(), key = lambda x: x[1][1], reverse=True)
    most_abund_seq_names = [tup[0] for tup in sorted_list]

    # from the above sorted list we can then plot these sequences with colour and all others with grey scale
    # lets make a coulour dictionary for the most common types
    colour_list = get_colour_list()
    colour_dict = {}
    num_coloured_seqs = 30

    # we will also need a grey palette for those sequences that are not the ones being annotated
    grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    # make the dict
    for i in range(len(most_abund_seq_names)):
        if i < num_coloured_seqs:
            colour_dict[most_abund_seq_names[i]] = colour_list[i]
        else:
            colour_dict[most_abund_seq_names[i]] = grey_palette[i % 6]

    # add the 'low' key and assign it to grey for later on
    # the low category will be created later on to hold the grouped abundances of sequences in samples
    # below a certain rel abund cutoff
    colour_dict['low'] = '#D0CFD4'

    # for plotting we will also need to collect the 'rectangle' objects from the plotting process to use to make
    # the legend which we will display on the plot right at the end
    # this is going to be a little tricky to collect as not every sample type group that we are plotting
    # is going to have all of the top n sequences. So, we will have to pick up rectangles when they come up in
    # the plotting.
    # to keep track of which sequences we need we have:
    top_n_seqs_to_get = most_abund_seq_names[:num_coloured_seqs]
    # to keep track of which sequences we have already collected we have:
    top_n_seqs_done = []
    #finally to collect the rectangles and store them in order despite collecting them out of order we have:
    legend_rectangle_holder = [[] for i in range(num_coloured_seqs)]


    # set up the plotting environment


    # rather than work with the subplot standard layout we will use subplot2grid for our set up and have some hard
    # coded axes

    this_fig = plt.figure(figsize=(10, 8))

    ax_list = [
        # coral
        (plt.subplot2grid((6, 67), (0, 0), colspan=67), 'coral'),

        # mucus
        (plt.subplot2grid((6, 67), (1, 0), colspan=66), 'mucus'),

        # sea water
        (plt.subplot2grid((6, 67), (2, 0 ),  colspan=20), 'sea_water'),

        # sed_close
        (plt.subplot2grid((6, 67), (3, 0 ), colspan=40), 'sed'),

        # turf,
        (plt.subplot2grid((6, 67), (4, 0 ), colspan= 37), 'turf'),

    ]



    # now we need to get the actual plotting information for each of the samples.
    # we can do this sample type by sample type
    # we will create a local dataframe that will be a sub set of the main sp output dataframe but will
    # only contain the samples of the given sample type. It will also eventually only contain the sequence
    # information for the sample type in question. Normally I would make a 2D list to hold the plotting
    # information but in this case having all of the information in a dataframe is working really well and
    # this is how I will do it in future.


    for ax_item in ax_list:
        try:
            env_sp_output_df_prop = pickle.load(open('env_sp_output_df_prop_{}_grouped_sed_grouped.pickle'.format(ax_item[1]), 'rb'))
            sorted_list_of_env_specific_seqs = pickle.load(open('sorted_list_of_env_specific_seqs_{}_grouped_sed_grouped.pickle'.format(ax_item[1]), 'rb'))
        except:
            sys.stdout.write('\nGenerating plotting info for samples {}\n'.format(ax_item[1]))
            # currently we have something like 4000 sequences to plot which is too may
            # I think it will be much easier if we group the sequences that are found below a certain
            # threshold. I think the best way to do this is to create a slice of the main df that
            # contain the information for the samples of the env_type only

            # get subset of the main dfs that contain only the env_type samples from the given date
            env_info_df = info_df[(info_df['environ type'] == ax_item[1])]

            # we need to drop any samples that have only zeros in their columns
            # https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame
            env_info_df = env_info_df[(env_info_df.T != 0).any()]


            env_sp_output_df = sp_output_df.loc[env_info_df.index.values.tolist()]

            # we need to drop any samples that have only zeros in their columns
            # https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame
            env_sp_output_df = env_sp_output_df[(env_sp_output_df.T != 0).any()]

            # append a 'low' columns to the env_sp_ouput_df and populate with 0s
            env_sp_output_df['low'] = 0
            # now make a proportions version of the df, rather than absolute counts
            env_sp_output_df_prop = env_sp_output_df[:].div(env_sp_output_df[:].sum(axis=1), axis=0)

            # now as we work our way through we can sum up the low sequences into this column
            # we can then check for 0 columns and drop these.

            # get a list of the sequences found in the collection of samples of the given type and order
            # them according to summed rel_abundance acorss all samples. This should be the order in which
            # the samples are plotted
            # at the same time we can get the info we need for plotting

            summed_seq_rel_abund_across_smpls_dict = {seq: 0 for seq in list(sp_output_df)}

            # we are also going to need to put the samples in an order that makes sense.
            # Ideally I would like to take the time to do a proper travelling salesman analysis
            # and I was thinking of putting in an ant colony system to do this but...
            # I'm not sure we've got time for that. So let's do something a little more simple which
            # should still be effective and look groovy
            # Lets sort according to the majority sequence. I.e. lets put the samples into groups that
            # are defined by what their majorty sequence is, then lets plot in order of those groups.
            # within the groups we will plot in the order of the most abundant rel abund first.
            # to get this information we will have a dict to collect it. The key of the dict will
            # be the majority sequence with the value being a list that contains tuples. One tuple
            # for each sample in the list which will contain the sample_name and the re abund of the maj seq

            # Dict for holding the sample_sorting info
            sample_sorting_info_dict = defaultdict(list)

            for sp_index in env_sp_output_df_prop.index.values.tolist():
                sys.stdout.write('\r{}'.format(sp_index))
                # we need to get the name of the most abundant sequence and its rel abund for each sample
                # for the sample_sorting_info_dict
                most_abund_seq_name = env_sp_output_df_prop.loc[sp_index].idxmax(axis='index')

                rel_abund_of_most_abund_seq = env_sp_output_df_prop.loc[sp_index, most_abund_seq_name]

                apples = 'asdf'
                sample_sorting_info_dict[most_abund_seq_name].append((sp_index, rel_abund_of_most_abund_seq))
                # Put its seq rel abundances ino the summed_seq_rel_abund_across... dict

                for non_zero_seq in env_sp_output_df_prop.loc[sp_index][env_sp_output_df_prop.loc[sp_index] > 0].index:  # not including the final 'low' columns (this will be zero)
                    val = env_sp_output_df_prop.loc[sp_index, non_zero_seq]
                    # be sure to count the value of this cell and using it in judging which are the most
                    # abundant sequences before we check whether to relegate it to the 'low' column
                    summed_seq_rel_abund_across_smpls_dict[non_zero_seq] += val
                    if val < 0.005:
                        env_sp_output_df_prop.loc[sp_index, 'low'] += val
                        env_sp_output_df_prop.loc[sp_index, non_zero_seq] = 0


            # here we can get a sorted sample list using the sample_sorting_info_dict
            sorted_sample_list = []
            # we want to work through the sample_sorting_info_dict by the longest lists first
            sorted_keys_for_sample_sort = \
                [tup[0] for tup in sorted(sample_sorting_info_dict.items(), key=lambda x: len(x[1]), reverse=True)]
            # for each of teh maj seq groups
            for sorted_key in sorted_keys_for_sample_sort:
                # now within each of these lists we want to order according to the rel_abundance of the sequences
                sorted_list_of_samples_of_group = [tup[0] for tup in sorted(sample_sorting_info_dict[sorted_key], key=lambda x: x[1], reverse=True)]
                sorted_sample_list.extend(sorted_list_of_samples_of_group)

            # now we should re-order the df so that it is in the sample order of sorted_sample_list
            env_sp_output_df_prop = env_sp_output_df_prop.reindex(sorted_sample_list)

            # here we have a dict that contains the abundances of the seqs for the coral samples
            # we will plot the coral samples' seqs in the order of the sequences in this dict
            sorted_list_of_env_specific_seqs_tup \
                = sorted(summed_seq_rel_abund_across_smpls_dict.items(), key=lambda x: x[1], reverse=True)
            sorted_list_of_env_specific_seqs = [tup[0] for tup in sorted_list_of_env_specific_seqs_tup]

            # Now we check for zero only columns and drop them
            # we also need to remove these columns from the sorted
            non_zero_cols = list(env_sp_output_df_prop.loc[:, (env_sp_output_df_prop != 0).any(axis=0)])
            sorted_list_of_env_specific_seqs = [seq for seq in sorted_list_of_env_specific_seqs if seq in non_zero_cols]

            # add the 'low' column only if the 'low' column is a non-zero column, i.e. used
            if 'low' in non_zero_cols:
                sorted_list_of_env_specific_seqs.append('low')

            # now drop the cols
            env_sp_output_df_prop = env_sp_output_df_prop[non_zero_cols]

            # we also have the plotting list which contains the info we will be plotting.

            # the plotting_list is currently organised in a different order to that of the sorted_list_of_env...
            # we need to change this order
            env_sp_output_df_prop = env_sp_output_df_prop[sorted_list_of_env_specific_seqs]

            # we now need to transpose this
            env_sp_output_df_prop = env_sp_output_df_prop.transpose()

            # pickle the info we've been plotting
            pickle.dump(env_sp_output_df_prop, open('env_sp_output_df_prop_{}_grouped_sed_grouped.pickle'.format(ax_item[1]), 'wb'))
            pickle.dump(sorted_list_of_env_specific_seqs, open('sorted_list_of_env_specific_seqs_{}_grouped_sed_grouped.pickle'.format(ax_item[1]), 'wb'))

        # here we finally have the plotting lists in the order of the sorted_list_of_env
        bottom = [0 for smp in list(env_sp_output_df_prop)]
        bar_list = []
        # for each sample
        plotting_indices = range(len(list(env_sp_output_df_prop)))

        # for each sequence
        list_of_seqs = env_sp_output_df_prop.index.values.tolist()
        for i in range(len(list_of_seqs)):

            sys.stdout.write('\r{}/{}'.format(i, len(list(env_sp_output_df_prop.iloc[:, 0]))))
            bar_list.append(ax_item[0].bar(plotting_indices, list(env_sp_output_df_prop.iloc[i]), 1, bottom,
                                                   color=colour_dict[sorted_list_of_env_specific_seqs[i]]))
            bottom = [L + M for L, M in zip(bottom, list(env_sp_output_df_prop.iloc[i]))]

            # check to see if the seq we are plotting is still in the top_n_seqs list
            # if it is in this list then we still need to grab a rectangle for plotting
            # the legend. Once we have grabbed the rectangle then we should remove the seq from the top_n_seqs list
            seq_name = list_of_seqs[i]

            if seq_name in top_n_seqs_to_get and seq_name not in top_n_seqs_done:
                # then this is a seq that we still need to get a rectangle for the legend
                legend_rectangle_holder[top_n_seqs_to_get.index(seq_name)].append(bar_list[-1][0])
                top_n_seqs_done.append(seq_name)

        # https://stackoverflow.com/questions/12998430/remove-xticks-in-a-matplotlib-plot
        ax_item[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        # https://stackoverflow.com/questions/15858192/how-to-set-xlim-and-ylim-for-a-subplot-in-matplotlib
        ax_item[0].set_xlim((-0.5, (len(list(env_sp_output_df_prop)) - 0.5)))
        # https://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib
        ax_item[0].spines['right'].set_visible(False)
        ax_item[0].spines['top'].set_visible(False)

        # To make the legend for this mo fo is going to be a little tricky. The legend argument basically takes
        # two lists. The first list should contain a rectangle object for each of the sequences (this
        # will be the coloured box). We get this object from the bar objects that we are creating.
        # The second list is a list of the labels. This should be easier. We can just use the
        # most_abund_seq_names[:30] for these.
        # to grab the rectangles, I think its best to pick them up during the plotting and hold them in a list
        # outside of the subplot loop. We will need a holder for the objects to populate.

        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticks.html#matplotlib.axes.Axes.set_xticks
        ax_item[0].set_yticks([1], minor=False)
        #https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html#matplotlib.axes.Axes.set_xticklabels
        ax_item[0].set_yticklabels(['1'])
        # if ax_item[2] == '22.08.2016':
        ax_item[0].set_ylabel(ax_item[1])


    #https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots/26892326
    plt.text(0, 0.5, 'ITS2 sequence relative abundance', va='center', rotation='vertical')
    ordered_rectange_list = [sub_list[0] for sub_list in legend_rectangle_holder]
    plt.legend(ordered_rectange_list, top_n_seqs_to_get, loc='lower center')
    #plt.tight_layout()
    plt.savefig('env_smp_diversity_bars_all_dates_grouped_sed_grouped.svg')
    plt.savefig('env_smp_diversity_bars_all_dates_grouped_sed_grouped.png')
    plt.show()

def create_quan_diversity_figs():
    ''' The purpose of these figs will be to compare the diversity found in the different sample types
    I envisage them being a sub plot each step of the QC:
     1 - Post QC
     2 - Symbiodinium
     3 - non-Symbiodinium
     4 - post - MED
     5 - pst MED to pre MED symbiodinium ratio
     number 5 will be important to consider and will hopefully hopefully be similar between all of the sample
     types. It is important because MED has the possibility to collapse diverstiy in a non-linear manner.
     e.g. in a sample with 1000 sequences, that has roughly the same spread as a sample with only 100 sequences
     it is possible that MED will collapse both of these samples into 10 meaningful sequences. However, we would
     then need to make sure to compare the different samples according to their pre-MED diversity. Although, we can
     also discuss when the post-MED diversity means.
     For each of teh sup plots I envisage there being a a 'set of plots' for each of the sample types
     For each of these sets would in trun contain a set of two plots. One for absolute and one for unique. Each one
     of these plots would be the actual datapoints on the left, and then a mean point with SD bars to the right of it
     . We can put the absoulte and unique values on different axes as the differences between these will be huge.
     We can achieve this by setting up subplot axes and then using the .errorbar and .scatter functions.
     In terms of collecting the data we should simply use the cleaned up dataframes that we created when making the
     diversity plots.'''

    sp_output_df = pickle.load(open('sp_output_df.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df.pickle', 'rb'))
    info_df = pickle.load(open('info_df.pickle', 'rb'))

    # lets make 5 subplot
    # according to the above categories

    f, axarr = plt.subplots(5, 1)
    # counter to reference which set of axes we are plotting on
    axarr_index = 0
    # y_axis_labels = ['raw_contigs', 'post_qc', 'Symbiodinium', 'non-Symbiodinium', 'post-MED', 'post-MED / pre-MED']
    y_axis_labels = ['raw_contigs',  'non-Symbiodinium','Symbiodinium',  'post-MED', 'post-MED / pre-MED']


    # cycle through these strings to help us with our conditionals
    # one of these for each of the subplots that we will create
    # we will make these useful tuples that will hold the actual name of the columns that the data we want will
    # be in so that we can pull these out of the dataframe easily
    # for sub_plot_type in [('raw_contigs',), ('post_qc_absolute_seqs', 'post_qc_unique_seqs'),
    #                       ('post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs'),
    #                       ('post_taxa_id_absolute_non_symbiodinium_seqs','post_taxa_id_unique_non_symbiodinium_seqs'),
    #                       ('post_med_absolute','post_med_unique'),
    #                       ('med_ratio', True) ]:
    for sub_plot_type in [('raw_contigs',),
                          ('post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs'),
                          ('post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs'),
                          ('post_med_absolute', 'post_med_unique'),
                          ('med_ratio', True)]:


        # for each of the sub plots we will want to grab the absolute and unique counts and plot these
        # for each of the sample types.
        # go environment type by environment type

        # we will create some x axis indicies to arranage where we will be ploting
        # we can be smart with these later on and create some nice spacing layouts but for the time
        # being lets just get things plotted. Let's have one idices for each sample type and work
        # relatively from there.
        ind = range(6)
        ind_index = 0


        if sub_plot_type[0] != 'raw_contigs': ax2 = axarr[axarr_index].twinx()

        axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
        env_types_list = ['coral', 'mucus', 'sea_water', 'sed_close', 'sed_far', 'turf']
        for env_type in env_types_list:


            if sub_plot_type[0] == 'raw_contigs':
                # here we will plot just the raw_contigs
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'turf':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')
                    axarr[axarr_index].spines['right'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    axarr[axarr_index].spines['top'].set_visible(False)

                    # set the ticks
                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    # set the xaxis title
                    axarr[axarr_index].set_xlabel('raw_contigs')

                    axarr[axarr_index].set_ylim((0, 1200000))

                ind_index += 1
            elif sub_plot_type[0] != 'med_ratio':
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')


                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = list(env_QC_info_df.loc[:, sub_plot_type[1]])
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    ax2.set_ylabel( '', color='r')
                    ax2.tick_params('y', colors='r')
                    axarr[axarr_index].spines['top'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    # ax2.spines['bottom'].set_visible(False)

                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
                    if sub_plot_type[0] == 'post_taxa_id_absolute_non_symbiodinium_seqs':
                        axarr[axarr_index].set_ylim((0, 150000))
                    elif sub_plot_type[0] == 'post_taxa_id_absolute_symbiodinium_seqs':
                        axarr[axarr_index].set_ylim((0, 800000))
                    elif sub_plot_type[0] == 'post_med_absolute':
                        axarr[axarr_index].set_ylim((0, 800000))

                ind_index += 1
            else:
                # here we need to ploto out the MED ratios.
                # these are simply going to be the med abosultes divided by the symbiodinium absolutes
                # and same for the uniques.
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples

                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_absolute']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_absolute_symbiodinium_seqs']))]
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')

                    axarr[axarr_index].spines['top'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    # ax2.spines['bottom'].set_visible(False)

                    axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].set_xlabel('Symbiodinium / post-MED')
                    axarr[axarr_index].set_ylim((0, 1.1))

                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_unique']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_unique_symbiodinium_seqs']))]
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    ax2.set_ylabel('', color='r')
                    ax2.tick_params('y', colors='r')

                ind_index += 1

        axarr_index += 1
    apples = 'asdf'
    f.text(0, 0.65, 'ITS2 absolute sequence abundance', va='center', rotation='vertical', color='b')
    f.text(1 - 0.01, 0.55, 'ITS2 unique sequence abundance', ha='center', va='center', rotation='vertical', color='r')
    f.text(0.07, 0.18, 'ratio', va='center', rotation='vertical', color='b')
    f.text(1 - 0.05, 0.18, 'ratio', va='center', rotation='vertical', color='r')


    plt.tight_layout()
    f.savefig('diversity_stats.svg')
    f.show()
    return

def create_quan_diversity_figs_all_dates():
    ''' This is a modification of the above so that all of the dates are considered.
    The purpose of these figs will be to compare the diversity found in the different sample types
    I envisage them being a sub plot each step of the QC:
     1 - Post QC
     2 - Symbiodinium
     3 - non-Symbiodinium
     4 - post - MED
     5 - pst MED to pre MED symbiodinium ratio
     number 5 will be important to consider and will hopefully hopefully be similar between all of the sample
     types. It is important because MED has the possibility to collapse diverstiy in a non-linear manner.
     e.g. in a sample with 1000 sequences, that has roughly the same spread as a sample with only 100 sequences
     it is possible that MED will collapse both of these samples into 10 meaningful sequences. However, we would
     then need to make sure to compare the different samples according to their pre-MED diversity. Although, we can
     also discuss when the post-MED diversity means.
     For each of teh sup plots I envisage there being a a 'set of plots' for each of the sample types
     For each of these sets would in trun contain a set of two plots. One for absolute and one for unique. Each one
     of these plots would be the actual datapoints on the left, and then a mean point with SD bars to the right of it
     . We can put the absoulte and unique values on different axes as the differences between these will be huge.
     We can achieve this by setting up subplot axes and then using the .errorbar and .scatter functions.
     In terms of collecting the data we should simply use the cleaned up dataframes that we created when making the
     diversity plots.'''

    sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
    info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))

    # remove sample P7-G07 as it has no Symbiodinium samples
    sp_output_df.drop('P7_G07', axis='index', inplace=True)
    QC_info_df.drop('P7_G07', axis='index', inplace=True)
    info_df.drop('P7_G07', axis='index', inplace=True)

    # lets make 5 subplot
    # according to the above categories

    f, axarr = plt.subplots(5, 1)
    # counter to reference which set of axes we are plotting on
    axarr_index = 0
    # y_axis_labels = ['raw_contigs', 'post_qc', 'Symbiodinium', 'non-Symbiodinium', 'post-MED', 'post-MED / pre-MED']
    y_axis_labels = ['raw_contigs',  'non-Symbiodinium','Symbiodinium',  'post-MED', 'post-MED / Symbiodinium']


    # cycle through these strings to help us with our conditionals
    # one of these for each of the subplots that we will create
    # we will make these useful tuples that will hold the actual name of the columns that the data we want will
    # be in so that we can pull these out of the dataframe easily
    for sub_plot_type in [('raw_contigs',),
                          ('post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs'),
                          ('post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs'),
                          ('post_med_absolute', 'post_med_unique'),
                          ('med_ratio', True)]:


        # for each of the sub plots we will want to grab the absolute and unique counts and plot these
        # for each of the sample types.
        # go environment type by environment type

        # we will create some x axis indicies to arranage where we will be ploting
        # we can be smart with these later on and create some nice spacing layouts but for the time
        # being lets just get things plotted. Let's have one idices for each sample type and work
        # relatively from there.
        ind = range(6)
        ind_index = 0


        if sub_plot_type[0] != 'raw_contigs': ax2 = axarr[axarr_index].twinx()

        axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
        env_types_list = ['coral', 'mucus', 'sea_water', 'sed_close', 'sed_far', 'turf']
        for env_type in env_types_list:


            if sub_plot_type[0] == 'raw_contigs':
                # here we will plot just the raw_contigs
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'turf':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')
                    axarr[axarr_index].spines['right'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    axarr[axarr_index].spines['top'].set_visible(False)

                    # set the ticks
                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    # set the xaxis title
                    axarr[axarr_index].set_xlabel('raw_contigs')

                    axarr[axarr_index].set_ylim((0, 1200000))

                ind_index += 1
            elif sub_plot_type[0] != 'med_ratio':
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')


                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = list(env_QC_info_df.loc[:, sub_plot_type[1]])
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    ax2.set_ylabel( '', color='r')
                    ax2.tick_params('y', colors='r')
                    axarr[axarr_index].spines['top'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    # ax2.spines['bottom'].set_visible(False)

                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
                    if sub_plot_type[0] == 'post_taxa_id_absolute_non_symbiodinium_seqs':
                        axarr[axarr_index].set_ylim((0, 150000))
                    elif sub_plot_type[0] == 'post_taxa_id_absolute_symbiodinium_seqs':
                        axarr[axarr_index].set_ylim((0, 800000))
                    elif sub_plot_type[0] == 'post_med_absolute':
                        axarr[axarr_index].set_ylim((0, 800000))

                ind_index += 1
            else:
                # here we need to ploto out the MED ratios.
                # these are simply going to be the med abosultes divided by the symbiodinium absolutes
                # and same for the uniques.
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples

                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_absolute']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_absolute_symbiodinium_seqs']))]
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')

                    axarr[axarr_index].spines['top'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    # ax2.spines['bottom'].set_visible(False)

                    axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].set_xlabel('post-MED / Symbiodinium')
                    axarr[axarr_index].set_ylim((0, 1.1))

                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_unique']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_unique_symbiodinium_seqs']))]
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    ax2.set_ylabel('', color='r')
                    ax2.tick_params('y', colors='r')

                ind_index += 1

        axarr_index += 1
    apples = 'asdf'
    f.text(0, 0.65, 'ITS2 absolute sequence abundance', va='center', rotation='vertical', color='b')
    f.text(1 - 0.01, 0.55, 'ITS2 unique sequence abundance', ha='center', va='center', rotation='vertical', color='r')
    f.text(0.07, 0.18, 'ratio', va='center', rotation='vertical', color='b')
    f.text(1 - 0.05, 0.18, 'ratio', va='center', rotation='vertical', color='r')


    plt.tight_layout()
    f.savefig('diversity_stats_all_dates.svg')
    f.show()
    return

#code for making the QC
def create_quan_diversity_figs_all_dates_sed_collapsed():
    ''' This is a modification of the above so that all of the dates are considered.
    The purpose of these figs will be to compare the diversity found in the different sample types
    I envisage them being a sub plot each step of the QC:
     1 - Post QC
     2 - Symbiodinium
     3 - non-Symbiodinium
     4 - post - MED
     5 - pst MED to pre MED symbiodinium ratio
     number 5 will be important to consider and will hopefully hopefully be similar between all of the sample
     types. It is important because MED has the possibility to collapse diverstiy in a non-linear manner.
     e.g. in a sample with 1000 sequences, that has roughly the same spread as a sample with only 100 sequences
     it is possible that MED will collapse both of these samples into 10 meaningful sequences. However, we would
     then need to make sure to compare the different samples according to their pre-MED diversity. Although, we can
     also discuss when the post-MED diversity means.
     For each of teh sup plots I envisage there being a a 'set of plots' for each of the sample types
     For each of these sets would in trun contain a set of two plots. One for absolute and one for unique. Each one
     of these plots would be the actual datapoints on the left, and then a mean point with SD bars to the right of it
     . We can put the absoulte and unique values on different axes as the differences between these will be huge.
     We can achieve this by setting up subplot axes and then using the .errorbar and .scatter functions.
     In terms of collecting the data we should simply use the cleaned up dataframes that we created when making the
     diversity plots.'''

    sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
    info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))

    # remove sample P7-G07 as it has no Symbiodinium samples
    sp_output_df.drop('P7_G07', axis='index', inplace=True)
    QC_info_df.drop('P7_G07', axis='index', inplace=True)
    info_df.drop('P7_G07', axis='index', inplace=True)

    # lets make 5 subplot
    # according to the above categories

    f, axarr = plt.subplots(5, 1, figsize=(10,8))
    # counter to reference which set of axes we are plotting on
    axarr_index = 0
    # y_axis_labels = ['raw_contigs', 'post_qc', 'Symbiodinium', 'non-Symbiodinium', 'post-MED', 'post-MED / pre-MED']
    y_axis_labels = ['raw_contigs',  'non-Symbiodinium','Symbiodinium',  'post-MED', 'post-MED / Symbiodinium']


    # cycle through these strings to help us with our conditionals
    # one of these for each of the subplots that we will create
    # we will make these useful tuples that will hold the actual name of the columns that the data we want will
    # be in so that we can pull these out of the dataframe easily
    for sub_plot_type in [('raw_contigs',),
                          ('post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs'),
                          ('post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs'),
                          ('post_med_absolute', 'post_med_unique'),
                          ('med_ratio', True)]:


        # for each of the sub plots we will want to grab the absolute and unique counts and plot these
        # for each of the sample types.
        # go environment type by environment type

        # we will create some x axis indicies to arranage where we will be ploting
        # we can be smart with these later on and create some nice spacing layouts but for the time
        # being lets just get things plotted. Let's have one idices for each sample type and work
        # relatively from there.
        ind = range(5)
        ind_index = 0


        if sub_plot_type[0] != 'raw_contigs':
            ax2 = axarr[axarr_index].twinx()

        axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])

        # we will convert the sed_close and sed_far to simply sed
        info_df = info_df.applymap(lambda x: 'sed' if "sed" in str(x) else x)
        env_types_list = ['coral', 'mucus', 'sea_water', 'sed', 'turf']
        for env_type in env_types_list:


            if sub_plot_type[0] == 'raw_contigs':
                # here we will plot just the raw_contigs
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'turf':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')
                    axarr[axarr_index].spines['right'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    axarr[axarr_index].spines['top'].set_visible(False)

                    # set the ticks
                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    # set the xaxis title
                    axarr[axarr_index].set_xlabel('raw_contigs')

                    axarr[axarr_index].set_ylim((0, 1200000))

                ind_index += 1
            elif sub_plot_type[0] != 'med_ratio':
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')


                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = list(env_QC_info_df.loc[:, sub_plot_type[1]])
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    ax2.set_ylabel( '', color='r')
                    ax2.tick_params('y', colors='r')
                    axarr[axarr_index].spines['top'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    # ax2.spines['bottom'].set_visible(False)

                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
                    if sub_plot_type[0] == 'post_taxa_id_absolute_non_symbiodinium_seqs':
                        axarr[axarr_index].set_ylim((0, 150000))
                    elif sub_plot_type[0] == 'post_taxa_id_absolute_symbiodinium_seqs':
                        axarr[axarr_index].set_ylim((0, 800000))
                    elif sub_plot_type[0] == 'post_med_absolute':
                        axarr[axarr_index].set_ylim((0, 800000))

                ind_index += 1
            else:
                # here we need to ploto out the MED ratios.
                # these are simply going to be the med abosultes divided by the symbiodinium absolutes
                # and same for the uniques.
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples

                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_absolute']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_absolute_symbiodinium_seqs']))]
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')

                    axarr[axarr_index].spines['top'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    # ax2.spines['bottom'].set_visible(False)

                    axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].set_xlabel('post-MED / Symbiodinium')
                    axarr[axarr_index].set_ylim((0, 1.1))

                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_unique']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_unique_symbiodinium_seqs']))]
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    ax2.set_ylabel('', color='r')
                    ax2.tick_params('y', colors='r')

                ind_index += 1

        axarr_index += 1
    apples = 'asdf'
    f.text(0, 0.65, 'ITS2 absolute sequence abundance', va='center', rotation='vertical', color='b')
    f.text(1 - 0.01, 0.55, 'ITS2 unique sequence abundance', ha='center', va='center', rotation='vertical', color='r')
    f.text(0.07, 0.18, 'ratio', va='center', rotation='vertical', color='b')
    f.text(1 - 0.05, 0.18, 'ratio', va='center', rotation='vertical', color='r')


    # plt.tight_layout()
    f.savefig('diversity_stats_all_dates_sed_grouped.svg')
    f.savefig('diversity_stats_all_dates_sed_grouped.png')
    f.show()
    return

def create_quan_diversity_figs_pre_MED_seqs():
    ''' This is a modification of the above so that all of the dates are considered.
    The purpose of these figs will be to compare the diversity found in the different sample types
    I envisage them being a sub plot each step of the QC:
     1 - Post QC
     2 - Symbiodinium
     3 - non-Symbiodinium
     4 - post - MED
     5 - pst MED to pre MED symbiodinium ratio
     number 5 will be important to consider and will hopefully hopefully be similar between all of the sample
     types. It is important because MED has the possibility to collapse diverstiy in a non-linear manner.
     e.g. in a sample with 1000 sequences, that has roughly the same spread as a sample with only 100 sequences
     it is possible that MED will collapse both of these samples into 10 meaningful sequences. However, we would
     then need to make sure to compare the different samples according to their pre-MED diversity. Although, we can
     also discuss when the post-MED diversity means.
     For each of teh sup plots I envisage there being a a 'set of plots' for each of the sample types
     For each of these sets would in trun contain a set of two plots. One for absolute and one for unique. Each one
     of these plots would be the actual datapoints on the left, and then a mean point with SD bars to the right of it
     . We can put the absoulte and unique values on different axes as the differences between these will be huge.
     We can achieve this by setting up subplot axes and then using the .errorbar and .scatter functions.
     In terms of collecting the data we should simply use the cleaned up dataframes that we created when making the
     diversity plots.'''

    sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
    info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))

    # remove sample P7-G07 as it has no Symbiodinium samples
    sp_output_df.drop('P7_G07', axis='index', inplace=True)
    QC_info_df.drop('P7_G07', axis='index', inplace=True)
    info_df.drop('P7_G07', axis='index', inplace=True)

    # lets make 5 subplot
    # according to the above categories

    f, axarr = plt.subplots(3, 1, figsize=(6,4))

    # # make the y axes logarithmic
    # for ax in axarr:
    #     ax.set_yscale('log')


    # counter to reference which set of axes we are plotting on
    axarr_index = 0
    # y_axis_labels = ['raw_contigs', 'post_qc', 'Symbiodinium', 'non-Symbiodinium', 'post-MED', 'post-MED / pre-MED']
    y_axis_labels = ['raw_contigs',  'non-Symbiodiniaceae','Symbiodiniaceae',  'post-MED', 'post-MED / Symbiodinium']


    # cycle through these strings to help us with our conditionals
    # one of these for each of the subplots that we will create
    # we will make these useful tuples that will hold the actual name of the columns that the data we want will
    # be in so that we can pull these out of the dataframe easily
    for sub_plot_type in [('raw_contigs',),
                          ('post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs'),
                          ('post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs')]:

        # The grids were confusing when there were two axes
        # axarr[axarr_index].grid(b=True, which='major', axis='y')
        # for each of the sub plots we will want to grab the absolute and unique counts and plot these
        # for each of the sample types.
        # go environment type by environment type

        # we will create some x axis indicies to arranage where we will be ploting
        # we can be smart with these later on and create some nice spacing layouts but for the time
        # being lets just get things plotted. Let's have one idices for each sample type and work
        # relatively from there.
        ind = range(5)
        ind_index = 0


        if sub_plot_type[0] != 'raw_contigs':
            ax2 = axarr[axarr_index].twinx()
            ax2.set_yscale('symlog')
            axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
            axarr[axarr_index].set_yscale('symlog')
        else:
            axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
            axarr[axarr_index].set_yscale('symlog')
            locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8), numticks=12)
            axarr[axarr_index].yaxis.set_minor_locator(locmin)
            axarr[axarr_index].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())




        # we will convert the sed_close and sed_far to simply sed
        info_df = info_df.applymap(lambda x: 'sed' if "sed" in str(x) else x)
        env_types_list = ['coral', 'mucus', 'sea_water', 'sed', 'turf']

        for env_type in env_types_list:


            if sub_plot_type[0] == 'raw_contigs':
                # here we will plot just the raw_contigs
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                # axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'turf':

                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')
                    axarr[axarr_index].spines['right'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    axarr[axarr_index].spines['top'].set_visible(False)

                    # set the ticks
                    # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    # axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    # set the xaxis title
                    axarr[axarr_index].set_xlabel('raw_contigs')


                    # axarr[axarr_index].set_ylim((0, 1200000))

                ind_index += 1
            elif sub_plot_type[0] != 'med_ratio':
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples
                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = list(env_QC_info_df.loc[:, sub_plot_type[0]])
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                # axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')


                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = list(env_QC_info_df.loc[:, sub_plot_type[1]])
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                # ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    if sub_plot_type[0] == 'post_taxa_id_absolute_symbiodinium_seqs':
                        ax2.set_ylabel('', color='r')
                        ax2.tick_params('y', colors='r')
                        axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                        axarr[axarr_index].set_xticklabels(env_types_list)
                        axarr[axarr_index].spines['top'].set_visible(False)
                        ax2.spines['top'].set_visible(False)
                        axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
                    else:
                        ax2.set_ylabel( '', color='r')
                        ax2.tick_params('y', colors='r')
                        axarr[axarr_index].spines['top'].set_visible(False)
                        # axarr[axarr_index].spines['bottom'].set_visible(False)
                        ax2.spines['top'].set_visible(False)
                        # ax2.spines['bottom'].set_visible(False)

                        # axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                        # axarr[axarr_index].set_xticklabels(env_types_list)
                        axarr[axarr_index].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                        ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                        axarr[axarr_index].set_xlabel(y_axis_labels[axarr_index])
                        # if sub_plot_type[0] == 'post_taxa_id_absolute_non_symbiodinium_seqs':
                        #     # axarr[axarr_index].set_ylim((0, 150000))
                        #     # axarr[axarr_index].set_yscale('symlog')
                        #     # formatter = FuncFormatter(log_10_product)
                        #     # axarr[0].xaxis.set_major_formatter(formatter)
                        #     # axarr[0].set_xlim(1e-1, 1e5)
                        #     # locmaj = matplotlib.ticker.LogLocator(base=10, numticks=12)
                        #     # axarr[axarr_index].yaxis.set_major_locator(locmaj)
                        #     # axarr[0].set_xticks = [10 ** i for i in range(2, 6)]
                        #     # axarr[0].set_xticklabels = [10 ** i for i in range(2, 6)]
                        #     # locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8), numticks=12)
                        #     # axarr[axarr_index].yaxis.set_minor_locator(locmin)
                        #     # axarr[axarr_index].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                        # elif sub_plot_type[0] == 'post_taxa_id_absolute_symbiodinium_seqs':
                        #     axarr[axarr_index].set_ylim((0, 800000))
                        # elif sub_plot_type[0] == 'post_med_absolute':
                        #     axarr[axarr_index].set_ylim((0, 800000))


                ind_index += 1
            else:
                # here we need to ploto out the MED ratios.
                # these are simply going to be the med abosultes divided by the symbiodinium absolutes
                # and same for the uniques.
                # get a sub df of the main df according to the env_type
                # get subset of the main dfs that contain only the coral samples

                env_info_df = info_df[info_df['environ type'] == env_type]
                env_QC_info_df = QC_info_df.loc[env_info_df.index.values.tolist()]
                sys.stdout.write('\nGenerating plotting info for {} samples in subplot type {}\n'
                                 .format(env_type, sub_plot_type))
                # the data we are going to be plotting is so simple that rather than collecting it and then
                # plotting it we may as well just go straight to plotting it from the df

                # PLOT ABSOLUTE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value
                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_absolute']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_absolute_symbiodinium_seqs']))]
                x_values = [ind[ind_index] for y in y_values]
                axarr[axarr_index].scatter(x_values, y_values, marker='.', s=1, c='b')

                # now plot the mean and error bars
                # I know there is a mean and SD function on a pandas series but it is throwing out all sorts of
                # erros so lest stick with what we know
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)
                axarr[axarr_index].scatter(x=ind[ind_index] + 0.125, y=mean, marker='s', s=8, c='b')
                # axarr[axarr_index].errorbar(x=ind[ind_index] + 0.125, y=mean, yerr=std, fmt='none', c='b')

                if env_type == 'coral':
                    axarr[axarr_index].set_ylabel('', color='b')
                    axarr[axarr_index].tick_params('y', colors='b')

                    axarr[axarr_index].spines['top'].set_visible(False)
                    # axarr[axarr_index].spines['bottom'].set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    # ax2.spines['bottom'].set_visible(False)

                    axarr[axarr_index].set_xticks([a + 0.1875 for a in range(6)], minor=False)
                    axarr[axarr_index].set_xticklabels(env_types_list)
                    axarr[axarr_index].set_xlabel('post-MED / Symbiodinium')
                    axarr[axarr_index].set_ylim((0, 1.1))

                # PLOT UNIQUE
                # first plot the actual datapoints
                # x will be the indices, y will be the actual value

                y_values = [tup[0] / tup[1] for tup in zip(list(env_QC_info_df.loc[:, 'post_med_unique']), list(
                    env_QC_info_df.loc[:, 'post_taxa_id_unique_symbiodinium_seqs']))]
                x_values = [ind[ind_index] + 0.250 for y in y_values]
                ax2.scatter(x_values, y_values, marker='.', s=1, c='r')

                # now plot the mean and error bars
                std = statistics.stdev(y_values)
                mean = statistics.mean(y_values)

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='o', s=8, c='r')
                # ax2.errorbar(x=ind[ind_index] + 0.375, y=mean, yerr=std, fmt='none', c='r')

                if env_type == 'coral':
                    ax2.set_ylabel('', color='r')
                    ax2.tick_params('y', colors='r')

                ind_index += 1

        axarr_index += 1
    apples = 'asdf'
    f.text(0.01, 0.55, 'absolute number of ITS2 sequences\n(log10)', va='center', ha='center', rotation='vertical', color='b')
    f.text(1 - 0.01, 0.40, 'unique number of ITS2 sequences\n(log10)', ha='center', va='center', rotation='vertical', color='r')
    # f.text(0.07, 0.18, 'ratio', va='center', rotation='vertical', color='b')
    # f.text(1 - 0.05, 0.18, 'ratio', va='center', rotation='vertical', color='r')


    plt.tight_layout()
    f.savefig('diversity_stats_no_MED.svg')
    f.savefig('diversity_stats_no_MED.png')
    f.show()
    return

# def make_venns():
#     ''' Here we will aim to make venn diagrams we will make a set of three 2 cirle venns which will
#     just show the comparison of the env_samples to the coral.
#     We will also make a set of four 3 circle venn diagrams which will show all 3 way combinations.
#     There is a neat little module called matplotlib_venn which we can use to do this and it is dead simple.
#     All you need to do is pass it a list of sets.'''
#
#     # read in the dataframes created previously
#     sp_output_df = pickle.load(open('sp_output_df.pickle', 'rb'))
#     QC_info_df = pickle.load(open('QC_info_df.pickle', 'rb'))
#     info_df = pickle.load(open('info_df.pickle', 'rb'))
#
#     # create a dictionary that is sample name to env_type
#     sample_name_to_sample_type_dict = {}
#     for info_index in info_df.index.values.tolist():
#         if info_df.loc[info_index, 'environ type'] == 'coral':
#             sample_name_to_sample_type_dict[info_index] = 'coral'
#         elif info_df.loc[info_index, 'environ type'] == 'sea_water':
#             sample_name_to_sample_type_dict[info_index] = 'sea_water'
#         elif info_df.loc[info_index, 'environ type'] == 'sed_far':
#             sample_name_to_sample_type_dict[info_index] = 'sed'
#         elif info_df.loc[info_index, 'environ type'] == 'sed_close':
#             sample_name_to_sample_type_dict[info_index] = 'sed'
#         elif info_df.loc[info_index, 'environ type'] == 'mucus':
#             sample_name_to_sample_type_dict[info_index] = 'mucus'
#         elif info_df.loc[info_index, 'environ type'] == 'turf':
#             sample_name_to_sample_type_dict[info_index] = 'turf'
#
#     # here we have the dict populated and we can now get to work pulling out the set of sequence
#     # names found in each of the sample types
#     sample_types = ['coral', 'sea_water', 'sed', 'mucus', 'turf']
#
#     set_list = [set() for type in sample_types]
#
#     # now work through the sp output and check to see which columns are non-zero columns and add these column
#     # labels into the respective sets
#     for sp_output_index in sp_output_df.index.values.tolist():
#         type_of_sample  = sample_name_to_sample_type_dict[sp_output_index]
#         non_zero_seqs = list(sp_output_df.loc[sp_output_index][sp_output_df.loc[sp_output_index] > 0].index)
#         set_list[sample_types.index(type_of_sample)].update(non_zero_seqs)
#
#     # here we should have the sets populated
#     # now create a dictionary that is the sample_type name to the set
#     sample_type_set_dict = {}
#     for smp_t in sample_types:
#         sample_type_set_dict[smp_t] = (smp_t, set_list[sample_types.index(smp_t)])
#
#     colour_dict = {
#         'coral' : 'orange',
#         'sed' : 'brown',
#         'sea_water' : 'blue',
#         'turf' : 'green',
#         'mucus' : 'gray'}
#
#
#     figure, axes = plt.subplots(2, 2)
#     ax_count = 0
#     # then simply do permutations for the 5 x 3
#
#     coral_combo_plot_list = [('coral', 'mucus'), ('coral', 'sea_water'), ('coral', 'turf'), ('coral', 'sed')]
#     for combo in coral_combo_plot_list:
#         set_list = [sample_type_set_dict[combo[0]][1], sample_type_set_dict[combo[1]][1]]
#         labels = [sample_type_set_dict[combo[0]][0], sample_type_set_dict[combo[1]][0]]
#         c = venn2(subsets=set_list, set_labels=labels, ax=axes[int(ax_count / 2)][ax_count % 2])
#         element_list = ['10', '01']
#         for i in range(2):
#             # for each path we need to work out which colour we want it to be
#             # we can do this with a simple dict outside of here
#             c.get_patch_by_id(element_list[i]).set_color(colour_dict[labels[i]])
#             c.get_patch_by_id(element_list[i]).set_edgecolor('none')
#             c.get_patch_by_id(element_list[i]).set_alpha(0.4)
#         ax_count += 1
#
#
#
#     figure.savefig('coral_combo_venn.svg')
#     figure.show()
#     plt.close()
#
#     figure, axes = plt.subplots(1, 1)
#     # now lets plot the other vann of the three env_types against each other
#     set_list = [sample_type_set_dict['sea_water'][1], sample_type_set_dict['turf'][1], sample_type_set_dict['sed'][1]]
#     labels = [sample_type_set_dict['sea_water'][0], sample_type_set_dict['turf'][0], sample_type_set_dict['sed'][0]]
#
#     c = venn3(subsets=set_list, set_labels=labels, ax=axes)
#     element_list = ['100', '010', '001']
#     for i in range(3):
#         # for each path we need to work out which colour we want it to be
#         # we can do this with a simple dict outside of here
#         c.get_patch_by_id(element_list[i]).set_color(colour_dict[labels[i]])
#         c.get_patch_by_id(element_list[i]).set_edgecolor('none')
#         c.get_patch_by_id(element_list[i]).set_alpha(0.4)
#
#     figure.savefig('env_combo_venn.svg')
#     figure.show()
#     apples = 'adf'
#
# def make_venns_all_dates():
#     ''' This is a modification of the above the only difference being that we are working with all samples
#     of all times. It would be nice to somehow represent what proportion of the total sequences are represented
#     by the unique samples but I don't actually think this is possible. Instead we will create labels that
#     show what proportion of the total sequences the unique sequences represent.
#     Here we will aim to make venn diagrams we will make a set of three 2 cirle venns which will
#     just show the comparison of the env_samples to the coral.
#     We will also make a set of four 3 circle venn diagrams which will show all 3 way combinations.
#     There is a neat little module called matplotlib_venn which we can use to do this and it is dead simple.
#     All you need to do is pass it a list of sets.'''
#
#     # read in the dataframes created previously
#     sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
#     QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
#     info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))
#
#
#     # in order to be able to calculate the proportion of a sample types sequencs that the various sequence regions
#     # represent in proportion to the total sequences we need to work with a prportion df
#     sp_output_df = sp_output_df[:].div(sp_output_df[:].sum(axis=1), axis=0)
#
#     # create a dictionary that is sample name to env_type
#     sample_name_to_sample_type_dict = {}
#     for info_index in info_df.index.values.tolist():
#         if info_df.loc[info_index, 'environ type'] == 'coral':
#             sample_name_to_sample_type_dict[info_index] = 'coral'
#         elif info_df.loc[info_index, 'environ type'] == 'sea_water':
#             sample_name_to_sample_type_dict[info_index] = 'sea_water'
#         elif info_df.loc[info_index, 'environ type'] == 'sed_far':
#             sample_name_to_sample_type_dict[info_index] = 'sed'
#         elif info_df.loc[info_index, 'environ type'] == 'sed_close':
#             sample_name_to_sample_type_dict[info_index] = 'sed'
#         elif info_df.loc[info_index, 'environ type'] == 'mucus':
#             sample_name_to_sample_type_dict[info_index] = 'mucus'
#         elif info_df.loc[info_index, 'environ type'] == 'turf':
#             sample_name_to_sample_type_dict[info_index] = 'turf'
#
#     # here we have the dict populated and we can now get to work pulling out the set of sequence
#     # names found in each of the sample types
#     sample_types = ['coral', 'sea_water', 'sed', 'mucus', 'turf']
#
#
#     set_list = [set() for type in sample_types]
#
#     # now work through the sp output and check to see which columns are non-zero columns and add these column
#     # labels into the respective sets
#     # here we will create a list of dictionaries, one per sample type and this will be
#     # key = sequence, value = cumulative rel_props
#     for sp_output_index in sp_output_df.index.values.tolist():
#         type_of_sample = sample_name_to_sample_type_dict[sp_output_index]
#         non_zero_seqs = list(sp_output_df.loc[sp_output_index][sp_output_df.loc[sp_output_index] > 0].index)
#         set_list[sample_types.index(type_of_sample)].update(non_zero_seqs)
#
#     # here we should have the sets populated
#     # now create a dictionary that is the sample_type name to the set
#     sample_type_set_dict = {}
#     for smp_t in sample_types:
#         sample_type_set_dict[smp_t] = (smp_t, set_list[sample_types.index(smp_t)])
#
#     colour_dict = {
#         'coral': 'orange',
#         'sed': 'brown',
#         'sea_water': 'blue',
#         'turf': 'green',
#         'mucus': 'gray'}
#
#     figure, axes = plt.subplots(2, 2)
#     ax_count = 0
#     # then simply do permutations for the 5 x 3
#
#     coral_combo_plot_list = [('coral', 'mucus'), ('coral', 'sea_water'), ('coral', 'turf'), ('coral', 'sed')]
#     for combo in coral_combo_plot_list:
#         set_list = [sample_type_set_dict[combo[0]][1], sample_type_set_dict[combo[1]][1]]
#         labels = [sample_type_set_dict[combo[0]][0], sample_type_set_dict[combo[1]][0]]
#         c = venn2(subsets=set_list, set_labels=labels, ax=axes[int(ax_count / 2)][ax_count % 2])
#         element_list = ['10', '01']
#         for i in range(2):
#             # for each path we need to work out which colour we want it to be
#             # we can do this with a simple dict outside of here
#             c.get_patch_by_id(element_list[i]).set_color(colour_dict[labels[i]])
#             c.get_patch_by_id(element_list[i]).set_edgecolor('none')
#             c.get_patch_by_id(element_list[i]).set_alpha(0.4)
#         ax_count += 1
#
#     figure.savefig('coral_combo_venn_all_dates.svg')
#     figure.show()
#     plt.close()
#
#
#
#
#     figure, axes = plt.subplots(1, 1)
#     # now lets plot the other vann of the three env_types against each other
#     set_list = [sample_type_set_dict['sea_water'][1], sample_type_set_dict['turf'][1],
#                 sample_type_set_dict['sed'][1]]
#     labels = [sample_type_set_dict['sea_water'][0], sample_type_set_dict['turf'][0], sample_type_set_dict['sed'][0]]
#
#     c = venn3(subsets=set_list, set_labels=labels, ax=axes)
#     element_list = ['100', '010', '001']
#     for i in range(3):
#         # for each path we need to work out which colour we want it to be
#         # we can do this with a simple dict outside of here
#         c.get_patch_by_id(element_list[i]).set_color(colour_dict[labels[i]])
#         c.get_patch_by_id(element_list[i]).set_edgecolor('none')
#         c.get_patch_by_id(element_list[i]).set_alpha(0.4)
#
#     figure.savefig('env_combo_venn_all_dates.svg')
#     figure.show()
#     apples = 'adf'

# code for making the venns
def make_venns_all_dates_rel_props_annotated():
    ''' This is a modification of the above the only difference being that we are working with all samples
    of all times. It would be nice to somehow represent what proportion of the total sequences are represented
    by the unique samples but I don't actually think this is possible. Instead we will create labels that
    show what proportion of the total sequences the unique sequences represent.
    Here we will aim to make venn diagrams we will make a set of three 2 cirle venns which will
    just show the comparison of the env_samples to the coral.
    We will also make a set of four 3 circle venn diagrams which will show all 3 way combinations.
    There is a neat little module called matplotlib_venn which we can use to do this and it is dead simple.
    All you need to do is pass it a list of sets.'''

    # read in the dataframes created previously
    sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
    info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))


    # in order to be able to calculate the proportion of a sample types sequencs that the various sequence regions
    # represent in proportion to the total sequences we need to work with a prportion df
    sp_output_df = sp_output_df[:].div(sp_output_df[:].sum(axis=1), axis=0)

    # create a dictionary that is sample name to env_type
    sample_name_to_sample_type_dict = {}
    for info_index in info_df.index.values.tolist():
        if info_df.loc[info_index, 'environ type'] == 'coral':
            sample_name_to_sample_type_dict[info_index] = 'coral'
        elif info_df.loc[info_index, 'environ type'] == 'sea_water':
            sample_name_to_sample_type_dict[info_index] = 'sea_water'
        elif info_df.loc[info_index, 'environ type'] == 'sed_far':
            sample_name_to_sample_type_dict[info_index] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'sed_close':
            sample_name_to_sample_type_dict[info_index] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'mucus':
            sample_name_to_sample_type_dict[info_index] = 'mucus'
        elif info_df.loc[info_index, 'environ type'] == 'turf':
            sample_name_to_sample_type_dict[info_index] = 'turf'

    # here we have the dict populated and we can now get to work pulling out the set of sequence
    # names found in each of the sample types
    sample_types = ['coral', 'sea_water', 'sed', 'mucus', 'turf']


    set_list = [set() for type in sample_types]

    # now work through the sp output and check to see which columns are non-zero columns and add these column
    # labels into the respective sets
    # here we will create a list of dictionaries, one per sample type and this will be
    # key = sequence, value = cumulative rel_props
    dict_of_info_dicts = {smp_type : defaultdict(float) for smp_type in sample_types}
    for sp_output_index in sp_output_df.index.values.tolist():
        sys.stdout.write('\r sample: {}'.format(sp_output_index))
        type_of_sample = sample_name_to_sample_type_dict[sp_output_index]

        non_zero_seqs = list(sp_output_df.loc[sp_output_index][sp_output_df.loc[sp_output_index] > 0].index)
        for seq_name in non_zero_seqs:
            dict_of_info_dicts[type_of_sample][seq_name] += sp_output_df.loc[sp_output_index, seq_name]
        set_list[sample_types.index(type_of_sample)].update(non_zero_seqs)

    # here we have both the set_lists and the dictionaries populated

    # here we should have the sets populated
    # now create a dictionary that is the sample_type name to the set
    sample_type_set_dict = {}
    for smp_t in sample_types:
        sample_type_set_dict[smp_t] = (smp_t, set_list[sample_types.index(smp_t)])

    colour_dict = {
        'coral': 'orange',
        'sed': 'brown',
        'sea_water': 'blue',
        'turf': 'green',
        'mucus': 'gray'}

    figure, axes = plt.subplots(2, 2, figsize=(14,10))
    ax_count = 0


    coral_combo_plot_list = [('coral', 'mucus'), ('coral', 'sea_water'), ('coral', 'turf'), ('coral', 'sed')]
    for combo in coral_combo_plot_list:
        set_list = [set([seq_name for seq_name in dict_of_info_dicts[combo[i]].keys()]) for i in range(2)]
        # set_list = [sample_type_set_dict[combo[0]][1], sample_type_set_dict[combo[1]][1]]
        labels = [combo[i] for i in range(2)]
        # labels = [sample_type_set_dict[combo[0]][0], sample_type_set_dict[combo[1]][0]]
        c = venn2(subsets=set_list, set_labels=labels, ax=axes[int(ax_count / 2)][ax_count % 2])
        v = venn2_circles(subsets=set_list, linestyle='dashed', color='grey', ax=axes[int(ax_count / 2)][ax_count % 2])
        element_list = ['10', '01']
        for i in range(2):
            # for each path we need to work out which colour we want it to be
            # we can do this with a simple dict outside of here
            c.get_patch_by_id(element_list[i]).set_color(colour_dict[labels[i]])
            c.get_patch_by_id(element_list[i]).set_edgecolor('none')
            c.get_patch_by_id(element_list[i]).set_alpha(0.4)

        element_list = ['10', '11', '01']
        elements_text_list = generate_elements_text_list([dict_of_info_dicts[combo[0]], dict_of_info_dicts[combo[1]]])
        for i in range(len(element_list)):
            c.get_label_by_id(element_list[i]).set_text('{}\n{}'.format(elements_text_list[i][0], elements_text_list[i][1]))

        ax_count += 1

    figure.savefig('coral_combo_venn_all_dates_labels.svg')
    figure.savefig('coral_combo_venn_all_dates_labels.png')
    figure.show()
    plt.close()




    figure, axes = plt.subplots(1, 1, figsize=(7,7))
    # now lets plot the other vann of the three env_types against each other

    set_list = [set([key for key in dict_of_info_dicts['sea_water'].keys()]), set([key for key in dict_of_info_dicts['turf'].keys()]),
                set([key for key in dict_of_info_dicts['sed'].keys()])]
    labels = ['sea_water', 'turf', 'sed']

    c = venn3(subsets=set_list, set_labels=labels, ax=axes)
    v = venn3_circles(subsets=set_list, linestyle='dashed', color='grey', ax=axes)
    element_list = ['100', '010', '001']
    for i in range(3):
        # for each path we need to work out which colour we want it to be
        # we can do this with a simple dict outside of here
        c.get_patch_by_id(element_list[i]).set_color(colour_dict[labels[i]])
        c.get_patch_by_id(element_list[i]).set_edgecolor('none')
        c.get_patch_by_id(element_list[i]).set_alpha(0.4)

    element_list = ['100', '110', '010', '011', '001', '101', '111']
    elements_text_list = generate_elements_text_list([dict_of_info_dicts['sea_water'], dict_of_info_dicts['turf'], dict_of_info_dicts['sed']])
    for i in range(len(element_list)):
        c.get_label_by_id(element_list[i]).set_text('{}\n{}'.format(elements_text_list[i][0], elements_text_list[i][1]))

    figure.savefig('env_combo_venn_all_dates_labels.svg')
    figure.savefig('env_combo_venn_all_dates_labels.png')
    figure.show()
    apples = 'adf'

def generate_elements_text_list(dict_list):
    # if the dict_list is only 2 elements long then we are working out a 2 venn
    # else we are working out a 3 venn

    if len(dict_list) == 2:
        dict_one_set = set(dict_list[0].keys())
        dict_two_set = set(dict_list[1].keys())
        dict_one_sum = sum(dict_list[0].values())
        dict_two_sum = sum(dict_list[1].values())
        union_list = \
            [
                dict_one_set.difference(dict_two_set),
                dict_one_set.intersection(dict_two_set),
                dict_two_set.difference(dict_one_set)
            ]
        list_of_all_seqs = list(set(dict_list[0].keys()) | set(dict_list[1].keys()))
        # return a list of labels in the order of 10 11 01
        label_list_info = [[0 for j in range(len(dict_list))] for i in range(3)]
        for seq_name in list_of_all_seqs:
            if seq_name in union_list[0]:
                # then it is only in 10
                label_list_info[0][0] += dict_list[0][seq_name]
            elif seq_name in union_list[1]:
                # then it is in both of the dicts (11)
                label_list_info[1][0] += dict_list[0][seq_name]
                label_list_info[1][1] += dict_list[1][seq_name]
            elif seq_name in union_list[2]:
                # then the seq is only in 01
                label_list_info[2][1] += dict_list[1][seq_name]

        # At this point we have the label list populated
        # now we can generate the label for eachof the segments
        # order: 10, 11, 01
        label_list_name = []
        # 10
        label_list_name.append((len(union_list[0]),'{:.2f}:{:.2f}'.format(label_list_info[0][0]/dict_one_sum, 0)))
        # 11
        label_list_name.append((len(union_list[1]), '{:.2f}:{:.2f}'.format(label_list_info[1][0]/dict_one_sum, label_list_info[1][1]/dict_two_sum)))
        # 01
        label_list_name.append((len(union_list[2]), '{:.2f}:{:.2f}'.format(0, label_list_info[2][1] / dict_two_sum)))
        return label_list_name
    else:
        ## here we are dealing with a venn 3
        # element_list = ['100', '110', '010', '011', '001', '101', '111']
        dict_one_set = set(dict_list[0].keys())
        dict_two_set = set(dict_list[1].keys())
        dict_three_set = set(dict_list[2].keys())
        dict_one_sum = sum(dict_list[0].values())
        dict_two_sum = sum(dict_list[1].values())
        dict_three_sum = sum(dict_list[2].values())

        union_list = \
            [
                #100
                dict_one_set.difference(dict_two_set).difference(dict_three_set),
                #110
                dict_one_set.intersection(dict_two_set).difference(dict_three_set),
                #010
                dict_two_set.difference(dict_one_set).difference(dict_three_set),
                #011
                dict_two_set.intersection(dict_three_set).difference(dict_one_set),
                #001
                dict_three_set.difference(dict_one_set).difference(dict_two_set),
                #101
                dict_one_set.intersection(dict_three_set).difference(dict_two_set),
                #111
                dict_one_set.intersection(dict_three_set).intersection(dict_two_set)
            ]
        list_of_all_seqs = list(set(dict_list[0].keys()) | set(dict_list[1].keys()) | set(dict_list[2].keys()))
        # return a list of labels in the order of '100', '110', '010', '011', '001', '101', '111'
        label_list_info = [[0 for j in range(len(dict_list))] for i in range(7)]
        for seq_name in list_of_all_seqs:
            if seq_name in union_list[0]:
                # then it is only in 100
                label_list_info[0][0] += dict_list[0][seq_name]
            elif seq_name in union_list[1]:
                # then it is in both of the dicts (110)
                label_list_info[1][0] += dict_list[0][seq_name]
                label_list_info[1][1] += dict_list[1][seq_name]
            elif seq_name in union_list[2]:
                # then the seq is only in 010
                label_list_info[2][1] += dict_list[1][seq_name]
            elif seq_name in union_list[3]:
                # then it is in both of the dicts (011)
                label_list_info[3][1] += dict_list[1][seq_name]
                label_list_info[3][2] += dict_list[2][seq_name]
            elif seq_name in union_list[4]:
                # then the seq is only in 001
                label_list_info[4][2] += dict_list[2][seq_name]
            elif seq_name in union_list[5]:
                # then it is in both of the dicts (101)
                label_list_info[5][0] += dict_list[0][seq_name]
                label_list_info[5][2] += dict_list[2][seq_name]
            elif seq_name in union_list[6]:
                # then the seq is in all three dicts 111
                label_list_info[6][0] += dict_list[0][seq_name]
                label_list_info[6][1] += dict_list[1][seq_name]
                label_list_info[6][2] += dict_list[2][seq_name]

        # At this point we have the label list populated
        # now we can generate the label for eachof the segments
        # order: '100', '110', '010', '011', '001', '101', '111'
        label_list_name = []
        # 100
        label_list_name.append((
                len(union_list[0]),
                '{:.2f}:{:.2f}:{:.2f}'.format(label_list_info[0][0] / dict_one_sum, 0, 0)
            ))
        # 110
        label_list_name.append((
            len(union_list[1]),
            '{:.2f}:{:.2f}:{:.2f}'.format(
                label_list_info[1][0] / dict_one_sum,
                label_list_info[1][1] / dict_two_sum,
                0)))
        # 010
        label_list_name.append((
            len(union_list[2]),
            '{:.2f}:{:.2f}:{:.2f}'.format(0, label_list_info[2][1] / dict_two_sum, 0)))
        # 011
        label_list_name.append((
            len(union_list[3]),
            '{:.2f}:{:.2f}:{:.2f}'.format(
                0,
                label_list_info[3][1] / dict_two_sum,
                label_list_info[3][2] / dict_three_sum)
        ))
        # 001
        label_list_name.append((
            len(union_list[4]),
            '{:.2f}:{:.2f}:{:.2f}'.format(0, 0, label_list_info[4][2] / dict_three_sum)))
        # 101
        label_list_name.append((
            len(union_list[5]),
            '{:.2f}:{:.2f}:{:.2f}'.format(
                label_list_info[5][0] / dict_one_sum,
                0,
                label_list_info[5][2] / dict_three_sum)
        ))
        # 111
        label_list_name.append((
            len(union_list[6]),
            '{:.2f}:{:.2f}:{:.2f}'.format(
                label_list_info[6][0] / dict_one_sum,
                label_list_info[6][1] / dict_two_sum,
                label_list_info[6][2] / dict_three_sum)))
        return label_list_name

def get_colour_list():
    colour_list = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5",
                  "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693",
                  "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
                  "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744",
                  "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68",
                  "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
                  "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
                  "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
                  "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
                  "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
                  "#C8D0F6", "#A3A489", "#806C66", "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59",
                  "#8ADBB4", "#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94",
                  "#7ED379", "#012C58", "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393",
                  "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                  "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5", "#E773CE",
                  "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4", "#00005F", "#A97399",
                  "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02",
                  "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
                  "#F4D749", "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9",
                  "#FFFFFE", "#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527",
                  "#8BB400", "#797868", "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C",
                  "#B88183", "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                  "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F", "#003109",
                  "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E", "#1A3A2A", "#494B5A",
                  "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700",
                  "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71", "#868E7E", "#98D058",
                  "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F",
                  "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
    return colour_list


# For making the rarefaction figure
def rarefaction_curves():
    ''' The purpose of this will be to create a figure that has a rarefaction curve for each of the sample types.
    I'm hoping that this will do an excellent job of mitigating the problem we currently have with regards to the
    differences in sequences returned from the different environments.
    We can look up a paper that supports the fact that relatively accurate inference can be made
    using the begining part of the curve.

    To draw the rarefaction curve we will work on a sample by sample basis. For each sample, we will create a list
    that contains the resundant seuqnces according to aboslute counts (massive list) and then simply do random
    selects on this without replacement to get an array of sequences in the pick order. We can then iterate through
    this picked order and put the seqs into a set. We can then get the len of the set at given interval. the interval
    vs. the set len will be our coordinates for the plotting. However, we will bootstrap this so we will do the above
    process per sample maybe 100 times. then for each of the intervals we can get an average. Then for all of the
    sampes of a given sample type we can average all of the diversity points at each interval and plot this point
    (maybe we also plot the individual points too). Obviously we will lose points as we work along our intervals.
    This is fine (and an inherent part of the data which we are precielcy trying to overcome with this rarefaction
    approach) but we should make the reader aware of this decrease by having an n on the graph. Also error bars would
    be nice. This may mean that we need to have different subplots for each of the curves though.
    '''

    # read in the dataframes
    sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
    info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))

    # remove sample P7-G07 as it has no Symbiodinium samples
    sp_output_df.drop('P7_G07', axis='index', inplace=True)
    QC_info_df.drop('P7_G07', axis='index', inplace=True)
    info_df.drop('P7_G07', axis='index', inplace=True)

    # remove samples that don't have any counts in the MED QC columns
    keeper_row_labels = QC_info_df.index[QC_info_df['post_med_absolute'] != 0]
    QC_info_df = QC_info_df.loc[keeper_row_labels]
    sp_output_df = sp_output_df.loc[keeper_row_labels]
    info_df = info_df.loc[keeper_row_labels]
    apples = 'asdf'

    # create a dictionary that is sample name to env_type
    sample_name_to_sample_type_dict = {}
    for info_index in info_df.index.values.tolist():
        if info_df.loc[info_index, 'environ type'] == 'coral':
            sample_name_to_sample_type_dict[info_index] = 'coral'
        elif info_df.loc[info_index, 'environ type'] == 'sea_water':
            sample_name_to_sample_type_dict[info_index] = 'sea_water'
        elif info_df.loc[info_index, 'environ type'] == 'sed_far':
            sample_name_to_sample_type_dict[info_index] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'sed_close':
            sample_name_to_sample_type_dict[info_index] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'mucus':
            sample_name_to_sample_type_dict[info_index] = 'mucus'
        elif info_df.loc[info_index, 'environ type'] == 'turf':
            sample_name_to_sample_type_dict[info_index] = 'turf'


    boot_iterations = 100

    # get a list of the sampling frequencies that we want to sample over
    sampling_frequencies = []
    additions_list = [0, 0.25, 0.5, 0.75]
    orders = range(1, 6)
    for order in orders:
        for addition in additions_list:
            sampling_frequencies.append(int(10 ** (order + addition)))

    if os.path.isfile('result_dict_rare_{}.pickle'.format(boot_iterations)):
        result_dict = pickle.load(open('result_dict_rare_{}.pickle'.format(boot_iterations), 'rb'))
        apples = 'asdf'
    else:

        # get a list of the sampling frequencies that we want to sample over
        sampling_frequencies = []
        additions_list = [0, 0.25, 0.5, 0.75]
        orders = range(1, 6)
        for order in orders:
            for addition in additions_list:
                sampling_frequencies.append(int(10 ** (order + addition)))

        # we will send one series to be bootstrapped to a core for MPing
        input_series_queue = Queue()

        num_proc = 26

        manager = Manager()
        result_dict = manager.dict()

        for smp_index in sp_output_df.index.values.tolist():
            sys.stdout.write('\rPrinting: {}'.format(smp_index))
            non_zero_indices = list(sp_output_df.loc[smp_index].nonzero()[0])
            non_zero_series = sp_output_df.loc[smp_index].iloc[non_zero_indices]

            with open('test_bob.txt', 'w') as f:
                for item in non_zero_series.index.values:
                    f.write('{}\n'.format(item))
            new_list = []
            with open('test_bob.txt', 'r') as f:
                new_list.extend([line.rstrip() for line in f])
            new_dict = {a:b for a, b, in zip(new_list, non_zero_series.values.tolist())}

            input_series_queue.put((smp_index, new_dict))

        for i in range(num_proc):
            input_series_queue.put('STOP')

        list_of_processes = []
        for N in range(num_proc):
            # p = Process(target=rarefaction_curve_worker, args=(input_series_queue, boot_iterations,
            #                                                    result_dict, sampling_frequencies))
            p = Process(target=rarefaction_curve_worker, args=(input_series_queue, boot_iterations,
                                                                result_dict, sampling_frequencies))
            list_of_processes.append(p)
            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(result_dict), open('result_dict_rare_{}.pickle'.format(boot_iterations), 'wb'))

    # we have our results that we can work with held in the result_dict
    # we should be able to work directly with this dictionary for plotting so lets set this up now

    colour_dict = {
        'coral': 'orange',
        'sed': 'brown',
        'sea_water': 'blue',
        'turf': 'green',
        'mucus': 'gray'}

    fig, axarr = plt.subplots(1, 6, sharey=True, sharex=True, figsize=(10,6))

    # axarr[0].set_xscale('log')

    axx_ind = 0
    env_type_list = ['coral', 'mucus', 'sea_water', 'sed', 'turf']

    # Dict that will hold a list of the averaged series grouped by the sample type
    data_info = defaultdict(list)

    # for each sample workout the means of the bottstraps as a series and add these to the dict
    for smp in result_dict.keys():
        env_type_of_smp = sample_name_to_sample_type_dict[smp]
        temp_df = pd.DataFrame(result_dict[smp])
        averaged_series = temp_df.mean(axis=0)
        averaged_series.name = smp
        data_info[env_type_of_smp].append(averaged_series)


    # here we should have all of the info we need to make a dataframe that can directly be used for plotting for
    # each of the environmental sample types

    # This will hold the pairs of x, y lists for each of the environment types so that we can plot them all on the
    # last plot
    line_holder = []
    for env_type in env_type_list:
        try:
            # env_type_df = pd.DataFrame(data_info[env_type], columns = sampling_frequencies)
            env_type_df = pd.DataFrame.from_items([(s.name, s) for s in data_info[env_type]]).T
        except:
            asdf = 'asdf'

        # here we need to add the columns to the df.
        # if we had points for all of the sampling frequencies then we can simply use the sampling frequency values
        # else we need to use a slice of the sampling_frequencies
        if len(list(env_type_df)) == len(sampling_frequencies):
            env_type_df.columns = sampling_frequencies
        else:
            env_type_df.columns = sampling_frequencies[:len(list(env_type_df))]

        mean_y = []
        mean_x = []
        num_samples = len(env_type_df.iloc[:,0])
        for col in list(env_type_df):
            # here plot the individual points (one pont for each of the samples of the env_type that have a point for
            # this sampling frequency
            y_list = env_type_df[col].dropna().tolist()
            x_list = [col for i in range(len(y_list))]
            axarr[axx_ind].scatter(x_list, y_list, marker='.' , color=colour_dict[env_type], s=1)
            axarr[axx_ind].text(x_list[0], max(y_list) + 10, str(len(y_list)), fontsize=8, ha='center')

            # only plot the line point if the number of samples remaining is > 1/3 of the total samples
            if len(y_list) < num_samples/3:
                break
            mean_y.append(sum(y_list)/len(y_list))
            mean_x.append(x_list[0])

        # now draw the line
        apples = 'asdf'
        axarr[axx_ind].plot(mean_x, mean_y, color=colour_dict[env_type])
        line_holder.append((mean_x, mean_y, colour_dict[env_type]))
        # plt.show()
        apples = 'asdf'

        axarr[axx_ind].set_xlabel(env_type)
        axx_ind += 1

    for tup in line_holder:
        axarr[5].plot(tup[0], tup[1], color=tup[2])
        axarr[5].set_xlabel('all')

    # ticks_list = [10 ** i for i in range(6)]
    # plt.set_xticks = [10 ** i for i in range(2, 6)]

    # NB getting the axes to behave well when logged was tricky but the below link gave a good solution

    axarr[0].set_xscale('log')
    # formatter = FuncFormatter(log_10_product)
    # axarr[0].xaxis.set_major_formatter(formatter)
    # axarr[0].set_xlim(1e-1, 1e5)
    locmaj = matplotlib.ticker.LogLocator(base=10, numticks=12)
    axarr[0].xaxis.set_major_locator(locmaj)
    # axarr[0].set_xticks = [10 ** i for i in range(2, 6)]
    # axarr[0].set_xticklabels = [10 ** i for i in range(2, 6)]
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8), numticks=12)
    axarr[0].xaxis.set_minor_locator(locmin)
    axarr[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    grid(True)
    # plt.xlabel('sampling point log10')
    # plt.ylabel('unique sequences')
    fig.text(0.5, 0.02, 'sampling point log10', ha='center', va='center')
    fig.text(0.06, 0.5, 'unique sequences', ha='center', va='center', rotation='vertical')
    # plt.tight_layout()
    plt.savefig('rarefaction_one_third_cutoff.svg')
    plt.savefig('rarefaction_one_third_cutoff.png')
    apples = 'asdf'

def rarefaction_curves_no_MED():
    '''I am going to modify this so that we do the same style of rarefaction only working with pre-MED sequences.'''

    ''' The purpose of this will be to create a figure that has a rarefaction curve for each of the sample types.
    I'm hoping that this will do an excellent job of mitigating the problem we currently have with regards to the
    differences in sequences returned from the different environments.
    We can look up a paper that supports the fact that relatively accurate inference can be made
    using the begining part of the curve.

    To draw the rarefaction curve we will work on a sample by sample basis. For each sample, we will create a list
    that contains the resundant seuqnces according to aboslute counts (massive list) and then simply do random
    selects on this without replacement to get an array of sequences in the pick order. We can then iterate through
    this picked order and put the seqs into a set. We can then get the len of the set at given interval. the interval
    vs. the set len will be our coordinates for the plotting. However, we will bootstrap this so we will do the above
    process per sample maybe 100 times. then for each of the intervals we can get an average. Then for all of the
    sampes of a given sample type we can average all of the diversity points at each interval and plot this point
    (maybe we also plot the individual points too). Obviously we will lose points as we work along our intervals.
    This is fine (and an inherent part of the data which we are precielcy trying to overcome with this rarefaction
    approach) but we should make the reader aware of this decrease by having an n on the graph. Also error bars would
    be nice. This may mean that we need to have different subplots for each of the curves though.
    '''

    # read in the dataframes
    sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
    info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))

    # remove sample P7-G07 as it has no Symbiodinium samples
    sp_output_df.drop('P7_G07', axis='index', inplace=True)
    QC_info_df.drop('P7_G07', axis='index', inplace=True)
    info_df.drop('P7_G07', axis='index', inplace=True)

    # remove samples that don't have any counts in the MED QC columns
    keeper_row_labels = QC_info_df.index[QC_info_df['post_taxa_id_absolute_symbiodinium_seqs'] != 0]
    QC_info_df = QC_info_df.loc[keeper_row_labels]
    sp_output_df = sp_output_df.loc[keeper_row_labels]
    info_df = info_df.loc[keeper_row_labels]
    apples = 'asdf'

    # create a dictionary that is sample name to env_type
    sample_name_to_sample_type_dict = {}
    for info_index in info_df.index.values.tolist():
        if info_df.loc[info_index, 'environ type'] == 'coral':
            sample_name_to_sample_type_dict[info_index] = 'coral'
        elif info_df.loc[info_index, 'environ type'] == 'sea_water':
            sample_name_to_sample_type_dict[info_index] = 'sea_water'
        elif info_df.loc[info_index, 'environ type'] == 'sed_far':
            sample_name_to_sample_type_dict[info_index] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'sed_close':
            sample_name_to_sample_type_dict[info_index] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'mucus':
            sample_name_to_sample_type_dict[info_index] = 'mucus'
        elif info_df.loc[info_index, 'environ type'] == 'turf':
            sample_name_to_sample_type_dict[info_index] = 'turf'

    boot_iterations = 100

    # get a list of the sampling frequencies that we want to sample over
    sampling_frequencies = []
    additions_list = [0, 0.25, 0.5, 0.75]
    orders = range(1, 6)
    for order in orders:
        for addition in additions_list:
            sampling_frequencies.append(int(10 ** (order + addition)))

    if os.path.isfile('result_dict_rarefaction_{}_no_MED.pickle'.format(boot_iterations)):
        result_dict = pickle.load(open('result_dict_rarefaction_{}_no_MED.pickle'.format(boot_iterations), 'rb'))
        apples = 'asdf'
    else:


        # we will send one series to be bootstrapped to a core for MPing
        input_series_queue = Queue()

        num_proc = 7

        manager = Manager()
        result_dict = manager.dict()

        # for each sample, create a tup that is the sample name and a dictionary
        # the dictionary should be the a sequence identifier as key and an absolute abundance as value
        # previously this was done using the sp_output_df. However now we will do this using the pre_MED_seq dump
        # once we have the sequence data the strucutre defined above we should be able to pass this into the
        # the rest of the script in order to do the actual modelling and plotting of the rarefaction curve

        list_of_sample_names_from_indices = sp_output_df.index.values.tolist()
        tuple_holding_list = create_smp_name_to_dict_tup_list(list_of_sample_names=list_of_sample_names_from_indices)

        # put each tup into the queue to process
        for smp_name_to_dict_tup in tuple_holding_list:
            input_series_queue.put(smp_name_to_dict_tup)

        for i in range(num_proc):
            input_series_queue.put('STOP')

        list_of_processes = []
        for N in range(num_proc):
            p = Process(target=rarefaction_curve_worker, args=(input_series_queue, boot_iterations,
                                                               result_dict, sampling_frequencies))
            list_of_processes.append(p)
            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(result_dict), open('result_dict_rarefaction_{}_no_MED.pickle'.format(boot_iterations), 'wb'))

    # we have our results that we can work with held in the result_dict
    # we should be able to work directly with this dictionary for plotting so lets set this up now

    colour_dict = {
        'coral': 'orange',
        'sed': 'brown',
        'sea_water': 'blue',
        'turf': 'green',
        'mucus': 'gray'}

    fig, axarr = plt.subplots(1, 6, sharey=True, sharex=True, figsize=(10, 6))

    # axarr[0].set_xscale('log')

    axx_ind = 0
    env_type_list = ['coral', 'mucus', 'sea_water', 'sed', 'turf']

    # Dict that will hold a list of the averaged series grouped by the sample type
    data_info = defaultdict(list)

    # for each sample workout the means of the bottstraps as a series and add these to the dict
    for smp in result_dict.keys():
        env_type_of_smp = sample_name_to_sample_type_dict[smp]
        temp_df = pd.DataFrame(result_dict[smp])
        averaged_series = temp_df.mean(axis=0)
        averaged_series.name = smp
        data_info[env_type_of_smp].append(averaged_series)

    # here we should have all of the info we need to make a dataframe that can directly be used for plotting for
    # each of the environmental sample types

    # This will hold the pairs of x, y lists for each of the environment types so that we can plot them all on the
    # last plot
    line_holder = []
    for env_type in env_type_list:
        try:
            # env_type_df = pd.DataFrame(data_info[env_type], columns = sampling_frequencies)
            env_type_df = pd.DataFrame.from_items([(s.name, s) for s in data_info[env_type]]).T
        except:
            asdf = 'asdf'

        # here we need to add the columns to the df.
        # if we had points for all of the sampling frequencies then we can simply use the sampling frequency values
        # else we need to use a slice of the sampling_frequencies
        if len(list(env_type_df)) == len(sampling_frequencies):
            env_type_df.columns = sampling_frequencies
        else:
            env_type_df.columns = sampling_frequencies[:len(list(env_type_df))]

        mean_y = []
        mean_x = []
        num_samples = len(env_type_df.iloc[:, 0])
        for col in list(env_type_df):
            # here plot the individual points (one pont for each of the samples of the env_type that have a point for
            # this sampling frequency
            y_list = env_type_df[col].dropna().tolist()
            x_list = [col for i in range(len(y_list))]
            axarr[axx_ind].scatter(x_list, y_list, marker='.', color=colour_dict[env_type], s=1)
            axarr[axx_ind].text(x_list[0], max(y_list) + 10, str(len(y_list)), fontsize=8, ha='center')

            # only plot the line point if the number of samples remaining is > 1/3 of the total samples
            if len(y_list) < num_samples / 3:
                break
            mean_y.append(sum(y_list) / len(y_list))
            mean_x.append(x_list[0])

        # now draw the line
        apples = 'asdf'
        axarr[axx_ind].plot(mean_x, mean_y, color=colour_dict[env_type])
        line_holder.append((mean_x, mean_y, colour_dict[env_type]))
        # plt.show()
        apples = 'asdf'

        axarr[axx_ind].set_xlabel(env_type)
        axx_ind += 1

    for tup in line_holder:
        axarr[5].plot(tup[0], tup[1], color=tup[2])
        axarr[5].set_xlabel('all')

    # ticks_list = [10 ** i for i in range(6)]
    # plt.set_xticks = [10 ** i for i in range(2, 6)]

    # NB getting the axes to behave well when logged was tricky but the below link gave a good solution

    axarr[0].set_xscale('log')
    # formatter = FuncFormatter(log_10_product)
    # axarr[0].xaxis.set_major_formatter(formatter)
    # axarr[0].set_xlim(1e-1, 1e5)
    locmaj = matplotlib.ticker.LogLocator(base=10, numticks=12)
    axarr[0].xaxis.set_major_locator(locmaj)
    # axarr[0].set_xticks = [10 ** i for i in range(2, 6)]
    # axarr[0].set_xticklabels = [10 ** i for i in range(2, 6)]
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8), numticks=12)
    axarr[0].xaxis.set_minor_locator(locmin)
    axarr[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    grid(True)
    # plt.xlabel('sampling point log10')
    # plt.ylabel('unique sequences')
    fig.text(0.5, 0.02, 'sampling point log10', ha='center', va='center')
    fig.text(0.06, 0.5, 'unique sequences', ha='center', va='center', rotation='vertical')
    # plt.tight_layout()
    plt.savefig('rarefaction_one_third_cutoff.svg')
    plt.savefig('rarefaction_one_third_cutoff.png')
    apples = 'asdf'


def create_smp_name_to_dict_tup_list(list_of_sample_names):
    pre_MED_seq_dump_dir = '/Users/humebc/Google Drive/projects/gabi_ITS2/pre_MED_seqs'
    sample_files = [f for f in os.listdir(pre_MED_seq_dump_dir)]
    tuple_holding_list = []
    for smp_f in sample_files:

        # find what the sample name is by parsing through the sample names from the df and finding which one fits in
        # we will use the count to make sure that only one name fits into the file and is therefore unique.
        count = 0
        for smp_index in list_of_sample_names:
            if smp_index in smp_f.replace('-', '_'):
                count += 1
                sample_name = smp_index

        if count > 1:
            sys.exit('Unique names not found')
        if count == 0:
            # then this may be one of the samples that we dropped and we can just continue
            continue
        print('Creating abundance dict for sample: {}'.format(sample_name))
        # now we have the sample name we need to create the dictionary. We can do this easily from the name file
        # we whould be sure to go through all of the clades and combine these into a single dictionary. At this
        # point we are just looking for different sequences and it doesn't matter what clade they are from

        # this is the directory that the different clade files will be in.
        # we will need to go through each of the clade directories for each of the samples and work with the
        # name file that we find in each one of these
        pre_MED_wkd = '{}/{}'.format(pre_MED_seq_dump_dir, smp_f)

        sample_clade_dirs = [f for f in os.listdir(pre_MED_wkd)]

        dict_to_populate = {}
        for clade_dir in sample_clade_dirs:
            pre_MED_wkd_clade = '{}/{}'.format(pre_MED_wkd, clade_dir)
            for smp_clade_file in os.listdir(pre_MED_wkd_clade):
                if '.name' in smp_clade_file:
                    # then this is the name file we want to be working with
                    nameFile = []
                    with open('{}/{}'.format(pre_MED_wkd_clade, smp_clade_file), 'r') as f:
                        nameFile.extend([line.rstrip() for line in f])

                    # here we have a clade's nameFile
                    # from this we can start to populate the dict_to_populate
                    for name_line in nameFile:
                        dict_to_populate[name_line.split('\t')[0]] = len(name_line.split('\t')[1].split(','))

        # At this point we should have been through each of the clade directories for a given sample
        # next we need to create the sample name to dictionary tuple and add this to the
        tuple_holding_list.append((sample_name, dict_to_populate))
    return tuple_holding_list


def rarefaction_curve_worker(input_queue, num_bootstraps, result_dict, sampling_frequencies):

    for name, working_dict in iter(input_queue.get, 'STOP'):

        sys.stdout.write('\n\nSample: {}'.format(name))
        # for each sample, perform the bootstrapping

        # the enormous list that will hold a non-redundant list of sequences for picking from
        picking_list = []
        for non_zero_col_labl in working_dict.keys():
            picking_list.extend([non_zero_col_labl for i in range(working_dict[non_zero_col_labl])])

        # this will hold the lists which will hold the results of a single bool
        # so this will hold the results of all of the bools
        sample_boot_result_holder = []
        for it_num in range(num_bootstraps):
            sys.stdout.write('\nbootstrap: {}\n'.format(it_num))
            # this is the set that we will populate to count number of unique seqs
            sample_result_holder = []
            temp_set = set()
            try:
                pick_array = np.random.choice(picking_list, len(picking_list), replace=False)
            except:
                asdf = 'asdf'
            for i in range(len(pick_array)):
                sys.stdout.write('\riteration: {}'.format(i))
                temp_set.add(pick_array[i])
                if i in sampling_frequencies:
                    # then we sample the length of the set
                    sample_result_holder.append(len(temp_set))
            sample_boot_result_holder.append(sample_result_holder)

        # here we have conducted the bootstrapping for one of the samples
        result_dict[name] = sample_boot_result_holder


rarefaction_curves_no_MED()
