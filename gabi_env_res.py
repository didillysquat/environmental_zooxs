import pandas as pd
import sys
import pickle
import os
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
import matplotlib
from collections import defaultdict, Counter
import sys
import statistics
import itertools
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import numpy as np
from multiprocessing import Queue, Process, Manager
from random import shuffle
from timeit import default_timer as timer
import random
import math
import subprocess
from plumbum import local
import networkx as nx

def convert_interleaved_to_sequencial_fasta_two(fasta_in):
    fasta_out = []
    for i in range(len(fasta_in)):
        if fasta_in[i].startswith('>'):
            if fasta_out:
                # if the fasta is not empty then this is not the first
                fasta_out.append(temp_seq_str)
            #else then this is the first sequence and there is no need to add the seq.
            temp_seq_str = ''
            fasta_out.append(fasta_in[i])
        else:
            temp_seq_str = temp_seq_str + fasta_in[i]
    #finally we need to add in the last sequence
    fasta_out.append(temp_seq_str)
    return fasta_out


def writeListToDestination(destination, listToWrite):
    # print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite) - 1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite) - 1:
                writer.write(listToWrite[i])
            i += 1

def readDefinedFileToList(filename):
    temp_list = []
    with open(filename, mode='r') as reader:
        temp_list = [line.rstrip() for line in reader]
    return temp_list

# CODE FOR PRODUCING THE STACKED BAR PLOTS
# TODO this will need producing with the no med data. We can maybe use the MED for the main publication and the
# no med can be in the supps.
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
    plt.savefig('stacked_bar_diversity.svg')
    plt.savefig('stacked_bar_diversity.png')
    plt.show()

# CODE FOR PRODUCING THE QC/DIVERSITY FIGURE
def create_quan_diversity_figs_no_MED_seqs():
    ''' This is the latest version of this figure where we are discounting anything to do with MED
    This will produce a figure with three subplots.
    Subplots will be for raw contigs, non-symbiodiniaceae and Symbiodiniacea sequences.
    Both the total number of sequences and the number of distinct sequence.'''

    sp_output_df = pickle.load(open('sp_output_df_all_dates.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df_all_dates.pickle', 'rb'))
    info_df = pickle.load(open('info_df_all_dates.pickle', 'rb'))

    # remove sample P7-G07 as it has no Symbiodinium samples
    sp_output_df.drop('P7_G07', axis='index', inplace=True)
    QC_info_df.drop('P7_G07', axis='index', inplace=True)
    info_df.drop('P7_G07', axis='index', inplace=True)

    # lets make 3 subplot
    # according to the above categories

    f, axarr = plt.subplots(3, 1, figsize=(6,4))


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


# CODE FOR PRODUCING THE TWO SETS OF VENNS (THIS WILL NEED MODIFICATION TO MAKE IT FOR THE NO MED SEQS)
# TODO this should also have a no_MED version of it made
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


# CODE FOR PRODUCING THE RAREFACTION CURVES ON A SAMPLE BY SAMPLE BASIS SPLIT BY ENV TYPE
#TODO add another subplot to the end of this which is a zoon in up to the  10^4 for comparisons sakes
# TODO create another figure that also include the modeled rarefaction curves
def rarefaction_curves_no_MED_sample():
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

    # I'm going to drop any of the samples that didn't have at least 100 sequences. Given that the 5th sampling
    # point is at 100 we drop any sample that has less than 5 sampling points

    result_dict = {k:v for k,v in result_dict.items() if len(v[0]) > 4}

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
            axarr[axx_ind].text(x_list[0], max(y_list) + 50, str(len(y_list)), fontsize=8, ha='center')

            # only plot the line point if the number of samples remaining is > 1/3 of the total samples
            if len(y_list) < 4:
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
    plt.savefig('rarefaction_sample_no_MED.svg')
    plt.savefig('rarefaction_sample_no_MED.png')
    apples = 'asdf'


# CODE FOR PRODUCING THE RAREFACTION CURVES ON AN ENV TYPE BY ENV TYPE BASIS
def rarefaction_curves_no_MED_community():
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
    sp_output_df = sp_output_df.loc[keeper_row_labels]
    info_df = info_df.loc[keeper_row_labels]

    env_type_list = ['coral', 'mucus', 'sea_water', 'sed', 'turf']

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
    orders = range(1, 8)
    for order in orders:
        for addition in additions_list:
            sampling_frequencies.append(int(10 ** (order + addition)))

    if os.path.isfile('result_dict_rarefaction_{}_no_MED_community.pickle'.format(boot_iterations)):
        result_dict = pickle.load(open('result_dict_rarefaction_{}_no_MED_community.pickle'.format(boot_iterations), 'rb'))
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
        tuple_holding_list = create_smp_name_to_dict_tup_list_community(list_of_sample_names=list_of_sample_names_from_indices)

        # here we can simply rework this tuple list to create a new tuple for each
        # of the environment types
        # To make this work we will need to change from working with seq names to the actual sequences
        # else we won't be able to compare between samples to see if the same sample has been counted in multiple
        # samples
        community_level_tup_list = []
        for index, env_type in enumerate(env_type_list):
            list_of_env_type_tuples = [tup for tup in tuple_holding_list if sample_name_to_sample_type_dict[tup[0]] == env_type]

            # here we have a list of tups that are of the env_type
            # https://stackoverflow.com/questions/11011756/is-there-any-pythonic-way-to-combine-two-dicts-adding-values-for-keys-that-appe/11011846
            counter_holder_env = Counter(dict())
            for tup in list_of_env_type_tuples:
                counter_holder_env += Counter(tup[1])

            community_level_tup_list.append((env_type, dict(counter_holder_env)))
            print('{}: len= {}'.format(env_type, sum([value for value in dict(counter_holder_env).values()])))

        # here we have the community_level_tup_list popualated
        # and we can now pass this into the below worker and it should still work in the same way

        # put each tup into the queue to process
        for env_type_to_dict_tup in community_level_tup_list:
            input_series_queue.put(env_type_to_dict_tup)

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

        pickle.dump(dict(result_dict), open('result_dict_rarefaction_{}_no_MED_community.pickle'.format(boot_iterations), 'wb'))

    # we have our results that we can work with held in the result_dict
    # we should be able to work directly with this dictionary for plotting so lets set this up now

    # I'm going to drop any of the samples that didn't have at least 100 sequences. Given that the 5th sampling
    # point is at 100 we drop any sample that has less than 5 sampling points

    result_dict = {k:v for k,v in result_dict.items() if len(v[0]) > 4}

    colour_dict = {
        'coral': 'orange',
        'sed': 'brown',
        'sea_water': 'blue',
        'turf': 'green',
        'mucus': 'gray'}

    fig, axarr = plt.subplots(1, 6, sharey=True, sharex=True, figsize=(10, 6))

    # axarr[0].set_xscale('log')

    axx_ind = 0


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
            axarr[axx_ind].text(x_list[0], max(y_list) + 50, str(len(y_list)), fontsize=8, ha='center')

            # only plot the line point if the number of samples remaining is > 1/3 of the total samples
            if len(y_list) < 4:
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
    plt.savefig('rarefaction_no_MED.svg')
    plt.savefig('rarefaction_no_MED.png')
    apples = 'asdf'


def create_smp_name_to_dict_tup_list(list_of_sample_names):
    pre_MED_seq_dump_dir = '/Users/humebc/Google_Drive/projects/gabi_ITS2/pre_MED_seqs'
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

def create_smp_name_to_dict_tup_list_community(list_of_sample_names, pre_MED_seq_dump_dir='/Users/humebc/Google_Drive/projects/gabi_ITS2/pre_MED_seqs'):

    sample_files = [f for f in os.listdir(pre_MED_seq_dump_dir)]
    tuple_holding_list = []
    for smp_f in sample_files:

        # find what the sample name is by parsing through the sample names from the df and finding which one fits in
        # we will use the count to make sure that only one name fits into the file and is therefore unique.
        count = 0
        for smp_index in list_of_sample_names:
            if smp_f.startswith(smp_index) or smp_f.replace('-', '_').startswith(smp_index):
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

                elif '.fasta' in smp_clade_file:
                    # then this is the name file we want to be working with
                    fastaFile = []
                    with open('{}/{}'.format(pre_MED_wkd_clade, smp_clade_file), 'r') as f:
                        fastaFile.extend([line.rstrip() for line in f])
                    # now check each of the fastas to make sure that we get rid of the extra 'A' if it is there
                    for i in range(1, len(fastaFile), 2):
                        if fastaFile[i].startswith('AATGG'):
                            fastaFile[i] = fastaFile[i][1:]
            # we need to make a fasta dict from the fasta
            fasta_dict = {fastaFile[i][1:]: fastaFile[i+1] for i in range(0, len(fastaFile),2)}

            # here we have a clade's nameFile
            # from this we can start to populate the dict_to_populate
            for name_line in nameFile:
                dict_to_populate[fasta_dict[name_line.split('\t')[0]]] = len(name_line.split('\t')[1].split(','))

        # At this point we should have been through each of the clade directories for a given sample
        # next we need to create the sample name to dictionary tuple and add this to the
        tuple_holding_list.append((sample_name, dict_to_populate))
    return tuple_holding_list


def rarefaction_curve_worker(input_queue, num_bootstraps, result_dict, sampling_frequencies):

    for name, working_dict in iter(input_queue.get, 'STOP'):
        random_time = 0
        shuffle_time = 0
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

            start = timer()
            pick_array = np.random.choice(picking_list, len(picking_list), replace=False)
            end = timer()
            random_time += end - start

            # start = timer()
            # shuffle(picking_list)
            # end = timer()
            # shuffle_time += end - start
            # sys.stdout.write('\n shuffle: {}, random: {}\n'.format(shuffle_time, random_time))

            for i in sampling_frequencies:
                if i < len(picking_list):
                    sys.stdout.write('\rsampling at: {}'.format(i))
                    sample_result_holder.append(len(set(picking_list[:i])))
            sample_boot_result_holder.append(sample_result_holder)

        # here we have conducted the bootstrapping for one of the samples
        result_dict[name] = sample_boot_result_holder


def extract_type_profiles():
    if os.path.isfile('/Users/humebc/Google_Drive/projects/gabi_ITS2/taxa_modelling//31_DBV_070918_2018-09-28_02-38-53.092387.profiles.absolute.txt'):
        type_abund_df = pd.read_csv('/Users/humebc/Google_Drive/projects/gabi_ITS2/taxa_modelling'
                                    '/31_DBV_070918_2018-09-28_02-38-53.092387.profiles.absolute.txt',
                                    sep='\t', lineterminator='\n', index_col=0)
    else:
        # then we are working on the remote server
        type_abund_df = pd.read_csv('{}/31_DBV_070918_2018-09-28_02-38-53.092387.profiles.absolute.txt'.format(os.getcwd()),
                                    sep='\t', lineterminator='\n', index_col=0)

    apples = 'asdf'

    # for each sample, (row) get the id of the maj type for the sample.
    # add this to a defaulat dict list where type is key and list of samples is list
    type_to_sample_list_dict = defaultdict(list)

    # create a smp_id to smp_name dict
    smp_id_to_name_dict = {index_id: type_abund_df.loc[index_id][0] for index_id in type_abund_df.index.values.tolist()[6:-16]}

    for smpl_id in type_abund_df.index.values.tolist()[6:-16]:
        # get the sample name series cutting off the first value which is the sample name. Convert type to int
        temp_series = type_abund_df.loc[smpl_id][1:].astype('int')
        # nonzero_indices = temp_series.nonzero()[0]
        # if len(nonzero_indices) == 1:
        # get id of max type
        max_type_id = temp_series.idxmax()
        # check that smp has a genuine max type e.g. not 0.
        max_value = temp_series[max_type_id]
        if max_value > 0:
            # add the sample to the type IDs list
            type_to_sample_list_dict[int(max_type_id)].append((smpl_id))

    # here we can work out which samples we should be looking for co-occuring seqs within for each
    # of the chosen types
    # I'm having to cut out the C3-C3gulf types (x2) because the ed sequencing is so low
    chosen_type_IDs = [10334, 10333, 9733, 9918, 10448, 9504, 10279, 9737, 10322, 10280, 9732, 10421, 10153, 10136]


    if os.path.isfile('type_profile_rel_abund_dict_holder_dict.pickle'):
        type_profile_rel_abund_dict_holder_dict = pickle.load(open('type_profile_rel_abund_dict_holder_dict.pickle', 'rb'))
    else:
        type_profile_rel_abund_dict_holder_dict = create_type_profile_rel_abund_dict(chosen_type_IDs,
                                                                                     smp_id_to_name_dict,
                                                                                     type_to_sample_list_dict)


def model_its_type_profiles():
    type_profile_rel_abund_dict_holder_dict = pickle.load(open('type_profile_rel_abund_dict_holder_dict.pickle', 'rb'))
    # Here we have the type_profile_rel_abund_dict_holder_dict populated with the profiles for each of our types
    # eventually we will want to graphically represent these using splitstree and networkx but for the time
    # being we want to do the modelling

    # first we can show that the seqs really do represent the defining seqs rath than just constant background.
    # if we compare the seqs found between two differnt types and see that the seqs are still high then
    # this could point towards background being maintained.
    # equally we should get rid of anysequences that are found in common to get a set of sequences that are truely
    # defining for each of the types. When we do the graphical representation we can highlight the seqs that we're really identifying

    # first lets get an idea of what kind of overlap we're looking at between the different types
    for type_tup_one, type_tup_two in itertools.combinations(list(type_profile_rel_abund_dict_holder_dict.items()), 2):
        # we simply want to see how many of the sequences are found in common, maybe keep an average
        print('{} and {} have {} seqs in common'.format(type_tup_one[0], type_tup_two[0], len(set(type_tup_one[1].keys()).intersection(set(type_tup_two[1].keys())))))

    type_to_clade_dict = {10334 : 'C',
                          10333: 'C',
                          9733: 'C',
                          9918: 'C',
                          10448: 'A',
                          9504: 'A',
                          10279: 'C',
                          9737: 'C',
                          10322: 'C',
                          10280: 'C',
                          9732: 'C',
                          10421: 'D',
                          10153: 'D',
                          10136: 'D'
                          }
    same_clade = []
    diff_clade = []

    # this is a default list dict that will be used to hold the sequences that need deleteing from a type
    # the sequences will be delted if they are found in common between two type profiles of the differnt clades
    # we will only delte the sequene from the type_profile that has it at the lower relative abundance as this is
    # likely the type profile and clade that it doesn't belong in.
    delete_wrong_clade_seqs(diff_clade, same_clade, type_profile_rel_abund_dict_holder_dict, type_to_clade_dict)

    # NB when we made the networks I thought that there might be an artefact as the most abundant seqs sometimes
    # weren't those with the radiations. However, I have checked this and it is not an artefact.
    draw_individual_type_networks(type_profile_rel_abund_dict_holder_dict)

    # here we want to draw the super network. This will be the network which is all of the types mixed together
    # just for diagramatic sake
    # we will make this network with a very high subsampling so that we can illustrate the loss of diversity
    # when we then do a much lower sampling, i.e. 1000.
    # we will somehow need to keep track of the proportion that each type contributes to every given sequence
    # this way we will be able to colour the super network so that we can see the mix in both the high
    # and low subsample level
    even_mix_type_dict = defaultdict(float)
    # we will have a dict dict here that will keep track of how much each type contributed to each sequence
    # key will be the sequence, and the value will be a list of tups which with each tup being
    # a type id and the relative proportion that the sequence was found in the type
    sequence_type_split_dict = defaultdict(list)
    for type_id_key, type_id_dict in type_profile_rel_abund_dict_holder_dict.items():
        for sequence_key, rel_abund_value in type_id_dict.keys():
            # be sure to add the sequence and the contrubution of this type to the sequence_type_split_dict too
            even_mix_type_dict[sequence_key] = rel_abund_value
            sequence_type_split_dict[sequence_key].append((type_id_key, rel_abund_value))


    # at this point we need to normalise the even mix dict
    total = sum(even_mix_type_dict.values())
    normalised_even_mix_type_dict = {k: v / total for k, v in even_mix_type_dict.items()}
    sys.stout.write('\rnew dict sums to {}'.format(sum(normalised_even_mix_type_dict.values())))
    # we should also normalise the sequence type split dict so that the proportions add up to one
    normalised_sequence_type_split_dict = {}
    for sequence_key, tup_list_value in sequence_type_split_dict.items():
        new_tup_list = []
        total = sum([tup[1] for tup in tup_list_value])
        new_tup_list = [(tup[0], tup[1]/total) for tup in tup_list_value]
        normalised_sequence_type_split_dict[sequence_key] = new_tup_list

    # now we should in theory be able to generate a splits tree network file with this

    splits_list_out_mix_1000 = generate_median_joining_network_from_dict_no_med(seq_abund_dict=normalised_even_mix_type_dict,
                                                                           type_id='even_mix', subsample_level=1000)

    splits_list_out_mix_100000 = generate_median_joining_network_from_dict_no_med(seq_abund_dict=normalised_even_mix_type_dict,
                                                                           type_id='even_mix', subsample_level=100000)
    # this will function slightly differently in that it will take into account that every node could have multiple
    # colours. If we can't use network x to do this for us then we can always do it by hand with circles.
    draw_network_split_colours()

    # here we have cleaned up type profiles
    # we have learnt that the high similarity between type profile are likely not due to random background zooxs
    # being carried through into type profiles as there is high similarity between similar types even from differnt studies
    # we also show very low similarity between type of different clade in the same studies which would again suggest that it is
    # a relation of types that gives the higher sequence similarity.
    # therefore we should believe the type profiles we have.

    # 1 - plot a rarefaction curve for each of the types on a single plot
    # also plot an eveness score as a product of sequencing depth
    # if we see that for the single types that the rarefaction curve looks like the turf then we can start to
    # say that the turf really is low. But if the rarefaction looks way higher then we know we have a problem
    # this could be a problem in itself.

    # plot_rarefaction_indi_gen_two(type_profile_rel_abund_dict_holder_dict)

    # 2 - for each type individually plot eveness as a function of sequencing depth
    plot_evenness(type_profile_rel_abund_dict_holder_dict)

    # 3 - for each type individually plot the average seq distance as a function of sequenceing depth

    # 4 then start the mixing: Do completely even and uneven where one predom, one quarter, rest small
    # do bootstraps on the even which represent just the sampling (i.e. one level of bootstrap)
    # but for the uneven need to change up which one is uneven, the within this do bootstrapping.
    # do this for eveness and average seq distance

    # this will show us whether we can use these as predictors.



    # 2 - calculate eveness
    # sampling frequecies


    # 3 - Caclulate average distance

    # maybe we can leave the distance stuff for a bit and get on to the modelling.

    # make a mix and calculate evenness
    # we can re use this code
    # instead of passing in the tuple for each of the type profiles that is id to dictionary
    # we can pass in tuples that are the mix id and a dictionary.
    # to get over biases associated with which types are in the mix, we should
    # pass in a random assortments of 1, 2, etc -->




    sampling_frequencies = []
    additions_list = [0, 0.25, 0.5, 0.75]
    orders = range(1, 6)
    for order in orders:
        for addition in additions_list:
            sampling_frequencies.append(int(10 ** (order + addition)))


    fig, axarr = plt.subplots(2, 4, figsize=(14, 6))

    # # GENERATE THE EVENESS DATA for the EVEN distribution
    list_of_mix_dictionaries_even = generate_mix_type_even(type_profile_rel_abund_dict_holder_dict)

    output_dict_shared_even_mix_evenness = dict(generate_evennessdata_for_even_distribution(list_of_mix_dictionaries_even, sampling_frequencies,
                                                                     type_profile_rel_abund_dict_holder_dict))

    # # GENERATE THE rarefaction data for the even distribution
    output_dict_shared_even_mix_rarefaction = dict(generate_rarefaction_data_mix_even(list_of_mix_dictionaries_even,
                                                                                      sampling_frequencies,
                                                                                      type_profile_rel_abund_dict_holder_dict))
    #
    # # PLOT Evenness
    #
    colour_dict = plot_evenness_even_dist(axarr, output_dict_shared_even_mix_evenness, sampling_frequencies)
    #

    #
    plot_rarefaction_mix_even(axarr, colour_dict, output_dict_shared_even_mix_rarefaction, sampling_frequencies)

    plot_rarefaction_mix_even_swap(axarr, output_dict_shared_even_mix_rarefaction, sampling_frequencies)
    #
    # # now we can move onto working with the uneven mixes. For this we will need to go back to
    # # we should be able to reuse a lot of code.

    list_of_mix_dictionaries_uneven = generate_mix_type_uneven(type_profile_rel_abund_dict_holder_dict)

    # here we have the dictionaries that represent the uneven mixes
    # we can then plot the same plots as before.
    # to do this we will need to do the eveness and rarefaction processing

    output_dict_shared_uneven_mix_rarefaction = dict(generate_rarefaction_data_mix_uneven(list_of_mix_dictionaries_uneven,
                                                                                 sampling_frequencies,
                                                                                 type_profile_rel_abund_dict_holder_dict))

    output_dict_shared_uneven_mix_evenness = dict(generate_evennessdata_for_uneven_distribution(list_of_mix_dictionaries_uneven, sampling_frequencies,
                                                                     type_profile_rel_abund_dict_holder_dict))

    colour_dict = plot_evenness_even_dist(axarr, output_dict_shared_uneven_mix_evenness, sampling_frequencies, uneven=True)

    plot_rarefaction_mix_even(axarr, colour_dict, output_dict_shared_uneven_mix_rarefaction, sampling_frequencies, uneven=True)

    plot_rarefaction_mix_even_swap(axarr, output_dict_shared_uneven_mix_rarefaction, sampling_frequencies, uneven=True)

    plt.tight_layout()
    apples = 'pears'


def draw_individual_type_networks(type_profile_rel_abund_dict_holder_dict):
    # we will pipe off here to test the network creation
    # to start with let's just work with trying to create a network from the first
    # dict
    # work out how many subplots we need
    num_types = len(type_profile_rel_abund_dict_holder_dict.items()) * 2
    for i in range(1, 10):
        if (i * i) > num_types:
            sub_plot_square = i
            break
    if os.path.isfile('network_colour_dict.pickle'):
        network_colour_dict = pickle.load(open('network_colour_dict.pickle', 'rb'))
    else:
        colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                              create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000,
                                                 num_cols=num_types,
                                                 time_out_iterations=10000)]
        network_colour_dict = {type_id: colour for type_id, colour in
                               zip(type_profile_rel_abund_dict_holder_dict.keys(), colour_palette_pas)}
        pickle.dump(network_colour_dict, open('network_colour_dict.pickle', 'wb'))

    f, axarr = plt.subplots(sub_plot_square, sub_plot_square, figsize=(10, 10))
    ax_count = 0
    for type_id, type_dict  in type_profile_rel_abund_dict_holder_dict.items():

        # get the ax object that we should be working with
        print('drawing network for ITS2 type profile {}'.format(type_id))


        # This returns two tuples. Each tuple represents a splitstree out file and an abundance dict
        # the first tup is currently the no med, and the second the med
        list_of_network_tups_to_draw = generate_median_joining_network_from_dict(type_dict, type_id, subsample_level=100000)
        for to_draw_tup in list_of_network_tups_to_draw:
            ax = axarr[int(ax_count / sub_plot_square)][ax_count % sub_plot_square]
            colour_for_network = network_colour_dict[type_id]
            draw_network(splits_tree_out_file=to_draw_tup[0], count_id_to_abund_dict=to_draw_tup[1], ax=ax, colour_for_network=colour_for_network, type_id=type_id)
            ax_count += 1
        apples = 'asdf'
    # remove the axes from the remaining two axes
    for i in range(2):
        ax = axarr[int(ax_count / sub_plot_square)][ax_count % sub_plot_square]
        ax.set_axis_off()
        ax_count += 1
    plt.savefig('individual_type_no_MED_networks_subsampled_with_med_M4_100000.svg')
    plt.savefig('individual_type_no_MED_networks_subsampled_with_med_M4_100000.png')


def delete_wrong_clade_seqs(diff_clade, same_clade, type_profile_rel_abund_dict_holder_dict, type_to_clade_dict):
    seqs_to_delete = defaultdict(list)
    for type_tup_one, type_tup_two in itertools.combinations(list(type_profile_rel_abund_dict_holder_dict.items()), 2):
        if type_to_clade_dict[type_tup_one[0]] == type_to_clade_dict[type_tup_two[0]]:
            same_clade.append(len(set(type_tup_one[1].keys()).intersection(set(type_tup_two[1].keys()))))
        else:
            diff_clade.append(len(set(type_tup_one[1].keys()).intersection(set(type_tup_two[1].keys()))))
            # get list of sequences that are found in common
            list_of_seqs_in_common_diff_clade = list(
                set(type_tup_one[1].keys()).intersection(set(type_tup_two[1].keys())))
            for seq in list_of_seqs_in_common_diff_clade:
                seq_rel_abund_one = type_tup_one[1][seq]
                seq_rel_abund_two = type_tup_two[1][seq]
                if seq_rel_abund_one > seq_rel_abund_two:
                    if seq not in seqs_to_delete[type_tup_two[0]]:
                        seqs_to_delete[type_tup_two[0]].append(seq)
                else:
                    if seq not in seqs_to_delete[type_tup_one[0]]:
                        seqs_to_delete[type_tup_one[0]].append(seq)
    print('av same_clade: {}'.format(sum(same_clade) / len(same_clade)))
    print('av diff_clade: {}'.format(sum(diff_clade) / len(diff_clade)))
    # delete the sequences that were found in common between the different clade comparisons.
    # delete the sequence from the type that had the lower abundance of it.
    for del_tup in seqs_to_delete.items():
        for seq_to_delete in del_tup[1]:
            type_id_in_Q = del_tup[0]
            del type_profile_rel_abund_dict_holder_dict[type_id_in_Q][seq_to_delete]


def create_type_profile_rel_abund_dict(chosen_type_IDs, smp_id_to_name_dict, type_to_sample_list_dict):
    pre_MED_seq_dir = '/Users/humebc/Google_Drive/projects/gabi_ITS2/taxa_modelling/pre_med_seqs'
    type_profile_rel_abund_dict_holder_dict = {}
    for t_id in chosen_type_IDs:
        # for each of the samples we will need to create a dict that is the sequence as key and the abundance of the sequence as value
        if t_id == 9737:
            apples = 'asdf'

        list_of_sample_names = [smp_id_to_name_dict[id] for id in type_to_sample_list_dict[t_id]]
        # if t_id == 9918:
        #     list_of_sample_names = [name for name in list_of_sample_names if name not in ['FS15SE7_FS15SE7_N704-S508', 'FS15SE8_FS15SE8_N705-S508']]
        if t_id == 9733:
            list_of_sample_names = [name for name in list_of_sample_names if name not in ['M_16_J0549_28-lost']]
        list_of_tups = create_smp_name_to_dict_tup_list_community(list_of_sample_names,
                                                                  pre_MED_seq_dump_dir=pre_MED_seq_dir)
        for tup in list_of_tups:
            print('{} has {} total sequences'.format(tup[0], sum(list(tup[1].values()))))
        # here we have a list of tups where key is sample name and value is a dict of seuence to abundance
        # for starters lets just see what happens if we just keep the sequences that are found in common
        # across all samples

        list_of_list_of_sequences = [list(tup[1].keys()) for tup in list_of_tups]

        common_seqs = set(list_of_list_of_sequences[0])

        print('Getting common seqs for type profile: {}'.format(t_id))
        for list_of_seqs in list_of_list_of_sequences[1:]:
            # we've got a little bit of a problem here in that some of the sequences have an extra A at the beginning
            # of them. They are causing us to find 0 seqs in common. However, it is quite binary
            # so we can just check and see if we have a total mis match and if so add an A to one of the parties
            print('Common seqs remaining: {}'.format(len(common_seqs)))
            common_seqs.intersection_update(list_of_seqs)
        print('Common seqs remaining: {}'.format(len(common_seqs)))
        common_seqs = list(common_seqs)
        # now that we have a list of the sequences that are found in common we should get the average relative
        # abundance of each of the sequences
        # Go back through each of the samples and get the relative abundance of each of the sequences
        # in the common list as a proportion of the other sequences in the common list
        # Then add these to a list of lists with each list containing the relative abundances for a given
        # sequence for a given sample
        relative_abundance_list_holder = [[] for seq in common_seqs]
        for tup in list_of_tups:
            smp_dict = tup[1]
            list_of_seq_absolute_abundances = []
            for common_seq in common_seqs:
                list_of_seq_absolute_abundances.append(smp_dict[common_seq])
            # here we have the list_of_seq_aboslute_abundances populated
            # now just work out these abundances as relative abundances and add to the relative_abundance_list_holder
            tot = sum(list_of_seq_absolute_abundances)
            for index, common_seq in enumerate(list_of_seq_absolute_abundances):
                relative_abundance_list_holder[index].append(common_seq / tot)

        # here we have the relative_abundance_list_holder populated.
        # we can now get an average of each of the lists and append these with their associated sequences names
        # to create a dictionary that is essentially our type profile for making the network from
        # as a sanity check we can also check that the sum of the dictionary values is equal to 1
        type_profile_dict = {}
        for index, rel_abund_list in enumerate(relative_abundance_list_holder):
            type_profile_dict[common_seqs[index]] = sum(rel_abund_list) / len(rel_abund_list)

        print('Relative abundance of type profile sum to: {}'.format(sum(type_profile_dict.values())))
        type_profile_rel_abund_dict_holder_dict[t_id] = type_profile_dict
    pickle.dump(type_profile_rel_abund_dict_holder_dict, open('type_profile_rel_abund_dict_holder_dict.pickle', 'wb'))
    return type_profile_rel_abund_dict_holder_dict


def generate_mix_type_even(type_profile_rel_abund_dict_holder_dict):
    if os.path.isfile('type_prof_mix_dictionary_list.pickle'):
        list_of_mix_dictionaries = pickle.load(open('type_prof_mix_dictionary_list.pickle', 'rb'))
    else:
        list_of_mix_dictionaries = []

        # first put in the singles
        for tup in type_profile_rel_abund_dict_holder_dict.items():
            list_of_mix_dictionaries.append(tup)

        # now do all of the others
        for i in range(2, len(type_profile_rel_abund_dict_holder_dict.items())):
            print('\nmaking the {} combo dicts'.format(i))
            counter = 0
            for combo in itertools.combinations(type_profile_rel_abund_dict_holder_dict.keys(), i):
                # create a counter that can hold the overall relative abundances for the combination of dictionaries
                combo_dict = Counter(dict())
                for single_dict_key in combo:
                    combo_dict += Counter(type_profile_rel_abund_dict_holder_dict[single_dict_key])
                # we then need to normalist the dict

                total = sum(combo_dict.values())
                normalised_dict = {k: v / total for k, v in combo_dict.items()}
                sys.stout.write('\rnew dict sums to {}'.format(sum(normalised_dict.values())))
                list_of_mix_dictionaries.append(('{}_{}'.format(i, counter), normalised_dict))
                counter += 1

        # here we have the list of dictionaries that we want to put into our MP and plot with
        pickle.dump(list_of_mix_dictionaries, open('type_prof_mix_dictionary_list.pickle', 'wb'))
    return list_of_mix_dictionaries

def generate_mix_type_uneven(type_profile_rel_abund_dict_holder_dict):
    # THis will produce us dicts that represent uneven mixes of Symbiodiniacea taxa
    # we will use a half life distribution so that relateive distributions are:
    # 1, 0.5, 0.25, 0.125, 0.0625, bottoming out at 0.0625
    distribution_list = [1, 0.5, 0.25, 0.125]
    list_of_low = [0.0625 for i in range(1000)]
    distribution_list.extend(list_of_low)

    if os.path.isfile('type_prof_mix_uneven_dictionary_list.pickle'):
        list_of_mix_dictionaries = pickle.load(open('type_prof_mix_uneven_dictionary_list.pickle', 'rb'))
    else:
        list_of_mix_dictionaries = []

        # first put in the singles
        # this will still be the same even for the unevens
        for tup in type_profile_rel_abund_dict_holder_dict.items():
            list_of_mix_dictionaries.append(tup)

        # now do all of the others
        for i in range(2, len(type_profile_rel_abund_dict_holder_dict.items())):
            print('\nmaking the {} combo dicts'.format(i))
            counter = 0
            for combo in itertools.combinations(type_profile_rel_abund_dict_holder_dict.keys(), i):
                # create a counter that can hold the overall relative abundances for the combination of dictionaries

                # firstly we should shuffle the combo list
                combo_shuffled = list(combo)
                random.shuffle(combo_shuffled)


                combo_dict = Counter(dict())
                for n, single_dict in enumerate(combo_shuffled):
                    # the dictionary abundances should be adjusted according to the distribution list
                    adjusted_dict = {k:v*distribution_list[n] for k, v in type_profile_rel_abund_dict_holder_dict[single_dict].items()}
                    combo_dict += Counter(adjusted_dict)
                # we then need to normalise the dict

                total = sum(combo_dict.values())
                normalised_dict = {k: v / total for k, v in combo_dict.items()}
                sys.stdout.write('\rnew dict sums to {}'.format(sum(normalised_dict.values())))
                list_of_mix_dictionaries.append(('{}_{}'.format(i, counter), normalised_dict))
                counter += 1

        # here we have the list of dictionaries that we want to put into our MP and plot with
        pickle.dump(list_of_mix_dictionaries, open('type_prof_mix_uneven_dictionary_list.pickle', 'wb'))
    return list_of_mix_dictionaries


def plot_rarefaction_mix_even(axarr, colour_dict, output_dict_shared_even_mix_rarefaction, sampling_frequencies, uneven=False):
    # list that will hold the series that will hold the averaged abundances at each bootstrap for each type profile
    series_holder = []
    # for each typeprofile workout the means of the bottstraps as a series and add these to the list
    for type_prof in output_dict_shared_even_mix_rarefaction.keys():
        temp_df = pd.DataFrame(output_dict_shared_even_mix_rarefaction[type_prof])
        averaged_series = temp_df.mean(axis=0)
        averaged_series.name = type_prof
        series_holder.append(averaged_series)
    # here we should have all of the info we need to make a dataframe that can directly be used for plotting
    # This will hold the pairs of x, y lists for each of the environment types so that we can plot them all on the
    # last plot
    env_type_df = pd.DataFrame.from_items([(s.name, s) for s in series_holder]).T
    # here we need to add the columns to the df.
    # use the sampling frequency values
    env_type_df.columns = sampling_frequencies
    # colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
    #                       create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000,
    #                                          num_cols=len(list(env_type_df)),
    #                                          time_out_iterations=10000)]
    # colour_dict = {col: colour for col, colour in zip(list(env_type_df), colour_palette_pas)}
    # plot a line for each sampling depth
    for col in list(env_type_df):
        x_val_line = []
        y_val_line = []
        for num_types in range(1, 14):
            # get a list of the keys that are for the number of samples num_types
            if num_types == 1:
                index_of = [k for k in env_type_df.index.values.tolist() if '_' not in str(k)]
            else:
                index_of = [k for k in env_type_df.index.values.tolist() if int(str(k).split('_')[0]) == num_types]

            # here plot the individual points (one pont for each of the samples of the env_type that have a point for
            # this sampling frequency
            y_list = [env_type_df.loc[ind][col] for ind in index_of]
            x_list = [num_types for i in y_list]
            # if uneven:
            #     axarr[1][1].scatter(x_list, y_list, color=colour_dict[col], marker='.')
            # else:
            #     axarr[0][1].scatter(x_list, y_list, color=colour_dict[col], marker='.')

            x_val_line.append(num_types)
            y_val_line.append(statistics.mean(y_list))
        # now draw the line
        if uneven:
            axarr[1][1].plot(x_val_line, y_val_line, color=colour_dict[col])

        else:
            axarr[0][1].plot(x_val_line, y_val_line, color=colour_dict[col])

    if uneven:
        axarr[1][1].set_yscale('log')
        axarr[1][1].set_ylabel('distinct sequences')
        axarr[1][1].set_xlabel('number of taxa')
    else:
        axarr[0][1].set_yscale('log')
        axarr[0][1].set_ylabel('distinct sequences')

def plot_rarefaction_mix_even_swap(axarr, output_dict_shared_even_mix_rarefaction, sampling_frequencies, uneven=False):
    # list that will hold the series that will hold the averaged abundances at each bootstrap for each type profile
    series_holder = []
    # for each typeprofile workout the means of the bottstraps as a series and add these to the list
    for type_prof in output_dict_shared_even_mix_rarefaction.keys():
        temp_df = pd.DataFrame(output_dict_shared_even_mix_rarefaction[type_prof])
        averaged_series = temp_df.mean(axis=0)
        averaged_series.name = type_prof
        series_holder.append(averaged_series)
    # here we should have all of the info we need to make a dataframe that can directly be used for plotting
    # This will hold the pairs of x, y lists for each of the environment types so that we can plot them all on the
    # last plot
    env_type_df = pd.DataFrame.from_items([(s.name, s) for s in series_holder]).T
    # here we need to add the columns to the df.
    # use the sampling frequency values
    env_type_df.columns = sampling_frequencies
    colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                          create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000,
                                             num_cols=len(range(1, 14)),
                                             time_out_iterations=10000)]
    colour_dict = {num_types: colour for num_types, colour in zip(range(1, 14), colour_palette_pas)}
    # plot a line for each sampling depth

    for num_types in range(1, 14):
        # get a list of the keys that are for the number of samples num_types
        if num_types == 1:
            index_of = [k for k in env_type_df.index.values.tolist() if '_' not in str(k)]
        else:
            index_of = [k for k in env_type_df.index.values.tolist() if int(str(k).split('_')[0]) == num_types]

        x_val_line = []
        y_val_line = []
        for col in list(env_type_df):

            # here plot the individual points (one pont for each of the samples of the env_type that have a point for
            # this sampling frequency
            y_list = [env_type_df.loc[ind][col] for ind in index_of]
            x_list = [col for i in y_list]
            # if uneven:
            #     axarr[1][2].scatter(x_list, y_list, color=colour_dict[num_types], marker='.')
            # else:
            #     axarr[0][2].scatter(x_list, y_list, color=colour_dict[num_types], marker='.')

            x_val_line.append(col)
            y_val_line.append(statistics.mean(y_list))
        # now draw the line
        if uneven:
            axarr[1][2].plot(x_val_line, y_val_line, color=colour_dict[num_types])
            axarr[1][2].set_xlabel('sequencing depth')
        else:
            axarr[0][2].plot(x_val_line, y_val_line, color=colour_dict[num_types])

    if uneven:
        axarr[1][2].set_xscale('log')
        axarr[1][2].set_ylabel('distinct sequences')
    else:
        axarr[0][2].set_xscale('log')
        axarr[0][2].set_ylabel('distinct sequences')

    for num_types in range(1, 14):
        # get a list of the keys that are for the number of samples num_types
        if num_types == 1:
            index_of = [k for k in env_type_df.index.values.tolist() if '_' not in str(k)]
        else:
            index_of = [k for k in env_type_df.index.values.tolist() if int(str(k).split('_')[0]) == num_types]

        x_val_line = []
        y_val_line = []
        for col in list(env_type_df)[:10]:

            # here plot the individual points (one pont for each of the samples of the env_type that have a point for
            # this sampling frequency
            y_list = [env_type_df.loc[ind][col] for ind in index_of]
            x_list = [col for i in y_list]
            # if uneven:
            #     axarr[1][3].scatter(x_list, y_list, color=colour_dict[num_types], marker='.')
            # else:
            #     axarr[0][3].scatter(x_list, y_list, color=colour_dict[num_types], marker='.')

            x_val_line.append(col)
            y_val_line.append(statistics.mean(y_list))
        # now draw the line
        if uneven:
            axarr[1][3].plot(x_val_line, y_val_line, color=colour_dict[num_types])

        else:
            axarr[0][3].plot(x_val_line, y_val_line, color=colour_dict[num_types])

    if uneven:
        axarr[1][3].set_xscale('log')
        axarr[1][3].set_ylabel('distinct sequences')
        axarr[1][3].set_xlabel('sequencing depth')
    else:
        axarr[0][3].set_xscale('log')
        axarr[0][3].set_ylabel('distinct sequences')
    apples  = 'asdf'

def generate_rarefaction_data_mix_even(list_of_mix_dictionaries, sampling_frequencies,
                                       type_profile_rel_abund_dict_holder_dict):
    if os.path.isfile('type_profile_rarefaction_mix_dict.pickle'):
        output_dict_shared_even_mix_rarefaction = pickle.load(open('type_profile_rarefaction_mix_dict.pickle', 'rb'))
    else:
        # create a dict to put the output in
        worker_manager = Manager()
        output_dict_shared_even_mix_rarefaction = worker_manager.dict()

        num_proc = 7

        # put in tups that are the items in the type_profile_rel_abund_dict_holder_dict
        input = Queue()

        list_of_processes = []
        # It will be too expensive to process all of the combinations, especially at this early stage.
        # for the time being lets just roll with 14 mixes for each number of samples
        # we should also run the same number of bootstraps on each to make the df easier to read.
        # we can always develop this more if it works out
        for i in range(1, (len(type_profile_rel_abund_dict_holder_dict.items()))):
            # get a list of the keys that are for the number of samples i
            if i == 1:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if '_' not in str(k[0])]
            else:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if int(str(k[0]).split('_')[0]) == i]

            # now we need to pick 14 of these keys randomly
            picking_list = list(
                np.random.choice(a=list_of_keys_of, size=len(type_profile_rel_abund_dict_holder_dict.items()), replace=False))

            # Here we have the list of keys to populate the worker with
            list_of_tups = [tup for tup in list_of_mix_dictionaries if tup[0] in picking_list]

            for tup in list_of_tups:
                input.put(tup)

        for i in range(num_proc):
            input.put('STOP')

        for N in range(num_proc):
            p = Process(target=rarefaction_curve_worker_indi_type_profile,
                        args=(input, 10, output_dict_shared_even_mix_rarefaction, sampling_frequencies))

            list_of_processes.append(p)

            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(output_dict_shared_even_mix_rarefaction),
                    open('type_profile_rarefaction_mix_dict.pickle', 'wb'))
    return output_dict_shared_even_mix_rarefaction

def generate_rarefaction_data_mix_uneven(list_of_mix_dictionaries, sampling_frequencies,
                                       type_profile_rel_abund_dict_holder_dict):
    if os.path.isfile('type_profile_rarefaction_uneven_mix_dict.pickle'):
        output_dict_shared_even_mix_rarefaction = pickle.load(open('type_profile_rarefaction_uneven_mix_dict.pickle', 'rb'))
    else:
        # create a dict to put the output in
        worker_manager = Manager()
        output_dict_shared_even_mix_rarefaction = worker_manager.dict()

        num_proc = 7

        # put in tups that are the items in the type_profile_rel_abund_dict_holder_dict
        input = Queue()

        list_of_processes = []
        # It will be too expensive to process all of the combinations, especially at this early stage.
        # for the time being lets just roll with 14 mixes for each number of samples
        # we should also run the same number of bootstraps on each to make the df easier to read.
        # we can always develop this more if it works out
        for i in range(1, (len(type_profile_rel_abund_dict_holder_dict.items()))):
            # get a list of the keys that are for the number of samples i
            if i == 1:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if '_' not in str(k[0])]
            else:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if int(str(k[0]).split('_')[0]) == i]

            # now we need to pick 14 of these keys randomly
            picking_list = list(
                np.random.choice(a=list_of_keys_of, size=len(type_profile_rel_abund_dict_holder_dict.items()), replace=False))

            # Here we have the list of keys to populate the worker with
            list_of_tups = [tup for tup in list_of_mix_dictionaries if tup[0] in picking_list]

            for tup in list_of_tups:
                input.put(tup)

        for i in range(num_proc):
            input.put('STOP')

        for N in range(num_proc):
            p = Process(target=rarefaction_curve_worker_indi_type_profile,
                        args=(input, 10, output_dict_shared_even_mix_rarefaction, sampling_frequencies))

            list_of_processes.append(p)

            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(output_dict_shared_even_mix_rarefaction),
                    open('type_profile_rarefaction_uneven_mix_dict.pickle', 'wb'))
    return output_dict_shared_even_mix_rarefaction


def plot_evenness_even_dist(axarr, output_dict_shared, sampling_frequencies, uneven=False):
    # list that will hold the series that will hold the averaged abundances at each bootstrap for each type profile
    series_holder = []
    # for each typeprofile workout the means of the bottstraps as a series and add these to the list
    for type_prof in output_dict_shared.keys():
        temp_df = pd.DataFrame(output_dict_shared[type_prof])
        averaged_series = temp_df.mean(axis=0)
        averaged_series.name = type_prof
        series_holder.append(averaged_series)
    # here we should have all of the info we need to make a dataframe that can directly be used for plotting
    # This will hold the pairs of x, y lists for each of the environment types so that we can plot them all on the
    # last plot
    env_type_df = pd.DataFrame.from_items([(s.name, s) for s in series_holder]).T
    # here we need to add the columns to the df.
    # use the sampling frequency values
    env_type_df.columns = sampling_frequencies
    colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                          create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000,
                                             num_cols=len(list(env_type_df)),
                                             time_out_iterations=10000)]
    colour_dict = {col: colour for col, colour in zip(list(env_type_df), colour_palette_pas)}
    # plot a line for each sampling depth
    for col in list(env_type_df):
        x_val_line = []
        y_val_line = []
        for num_types in range(1, 14):
            # get a list of the keys that are for the number of samples num_types
            if num_types == 1:
                index_of = [k for k in env_type_df.index.values.tolist() if '_' not in str(k)]
            else:
                index_of = [k for k in env_type_df.index.values.tolist() if int(str(k).split('_')[0]) == num_types]

            # here plot the individual points (one pont for each of the samples of the env_type that have a point for
            # this sampling frequency
            y_list = [env_type_df.loc[ind][col] for ind in index_of]
            x_list = [num_types for i in y_list]
            # if uneven:
            #     axarr[1][0].scatter(x_list, y_list, color=colour_dict[col], marker='.')
            # else:
            #     axarr[0][0].scatter(x_list, y_list, color=colour_dict[col], marker='.')

            x_val_line.append(num_types)
            y_val_line.append(statistics.mean(y_list))
        # now draw the line
        if uneven:
            axarr[1][0].plot(x_val_line, y_val_line, color=colour_dict[col])
            axarr[1][0].set_xlabel('number of taxa')
            axarr[1][0].set_ylabel('corrected Shannon equality')
        else:
            axarr[0][0].plot(x_val_line, y_val_line, color=colour_dict[col])
            axarr[0][0].set_ylabel('corrected Shannon equality')
    # axarr[0].set_xscale('log')
    return colour_dict




def generate_evennessdata_for_even_distribution(list_of_mix_dictionaries, sampling_frequencies,
                                                type_profile_rel_abund_dict_holder_dict):
    if os.path.isfile('type_profile_evenness_mix_dict.pickle'):
        output_dict_shared = pickle.load(open('type_profile_evenness_mix_dict.pickle', 'rb'))
    else:
        # create a dict to put the output in
        worker_manager = Manager()
        output_dict_shared = worker_manager.dict()

        num_proc = 7

        # put in tups that are the items in the type_profile_rel_abund_dict_holder_dict
        input = Queue()

        list_of_processes = []
        # It will be too expensive to process all of the combinations, especially at this early stage.
        # for the time being lets just roll with 14 mixes for each number of samples
        # we should also run the same number of bootstraps on each to make the df easier to read.
        # we can always develop this more if it works out
        for i in range(1, (len(type_profile_rel_abund_dict_holder_dict.items()))):
            # get a list of the keys that are for the number of samples i
            if i == 1:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if '_' not in str(k[0])]
            else:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if int(str(k[0]).split('_')[0]) == i]

            # now we need to pick 14 of these keys randomly
            picking_list = list(
                np.random.choice(a=list_of_keys_of, size=len(type_profile_rel_abund_dict_holder_dict.items()), replace=False))

            # Here we have the list of keys to populate the worker with
            list_of_tups = [tup for tup in list_of_mix_dictionaries if tup[0] in picking_list]

            for tup in list_of_tups:
                input.put(tup)

        for i in range(num_proc):
            input.put('STOP')

        for N in range(num_proc):
            p = Process(target=evenness_worker_indi_type_profile,
                        args=(input, 10, output_dict_shared, sampling_frequencies))

            list_of_processes.append(p)

            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(output_dict_shared), open('type_profile_evenness_mix_dict.pickle', 'wb'))
    return output_dict_shared

def generate_evennessdata_for_uneven_distribution(list_of_mix_dictionaries, sampling_frequencies,
                                                type_profile_rel_abund_dict_holder_dict):
    if os.path.isfile('type_profile_evenness_mix_uneven_dict.pickle'):
        output_dict_shared = pickle.load(open('type_profile_evenness_mix_uneven_dict.pickle', 'rb'))
    else:
        # create a dict to put the output in
        worker_manager = Manager()
        output_dict_shared = worker_manager.dict()

        num_proc = 7

        # put in tups that are the items in the type_profile_rel_abund_dict_holder_dict
        input = Queue()

        list_of_processes = []
        # It will be too expensive to process all of the combinations, especially at this early stage.
        # for the time being lets just roll with 14 mixes for each number of samples
        # we should also run the same number of bootstraps on each to make the df easier to read.
        # we can always develop this more if it works out
        for i in range(1, (len(type_profile_rel_abund_dict_holder_dict.items()))):
            # get a list of the keys that are for the number of samples i
            if i == 1:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if '_' not in str(k[0])]
            else:
                list_of_keys_of = [k[0] for k in list_of_mix_dictionaries if int(str(k[0]).split('_')[0]) == i]

            # now we need to pick 14 of these keys randomly
            picking_list = list(
                np.random.choice(a=list_of_keys_of, size=len(type_profile_rel_abund_dict_holder_dict.items()), replace=False))

            # Here we have the list of keys to populate the worker with
            list_of_tups = [tup for tup in list_of_mix_dictionaries if tup[0] in picking_list]

            for tup in list_of_tups:
                input.put(tup)

        for i in range(num_proc):
            input.put('STOP')

        for N in range(num_proc):
            p = Process(target=evenness_worker_indi_type_profile,
                        args=(input, 10, output_dict_shared, sampling_frequencies))

            list_of_processes.append(p)

            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(output_dict_shared), open('type_profile_evenness_mix_uneven_dict.pickle', 'wb'))
    return output_dict_shared


def plot_evenness(type_profile_rel_abund_dict_holder_dict):
    sampling_frequencies = []
    additions_list = [0, 0.25, 0.5, 0.75]
    orders = range(1, 6)
    for order in orders:
        for addition in additions_list:
            sampling_frequencies.append(int(10 ** (order + addition)))
    if os.path.isfile('type_profile_evenness_dict.pickle'):
        output_dict_shared = pickle.load(open('type_profile_evenness_dict.pickle', 'rb'))
    else:
        # create a dict to put the output in
        worker_manager = Manager()
        output_dict_shared = worker_manager.dict()

        num_proc = 7

        # put in tups that are the items in the type_profile_rel_abund_dict_holder_dict
        input = Queue()
        for tup in type_profile_rel_abund_dict_holder_dict.items():
            input.put(tup)

        for i in range(num_proc):
            input.put('STOP')

        list_of_processes = []
        for N in range(num_proc):
            p = Process(target=evenness_worker_indi_type_profile,
                        args=(input, 100, output_dict_shared, sampling_frequencies))
            list_of_processes.append(p)
            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(output_dict_shared), open('type_profile_evenness_dict.pickle', 'wb'))
    # PLOT Evenness
    fig, axarr = plt.subplots(1, 3, figsize=(10, 6))
    # list that will hold the series that will hold the averaged abundances at each bootstrap for each type profile
    series_holder = []
    # for each typeprofile workout the means of the bottstraps as a series and add these to the list
    for type_prof in output_dict_shared.keys():
        temp_df = pd.DataFrame(output_dict_shared[type_prof])
        averaged_series = temp_df.mean(axis=0)
        averaged_series.name = type_prof
        series_holder.append(averaged_series)
    # here we should have all of the info we need to make a dataframe that can directly be used for plotting
    # This will hold the pairs of x, y lists for each of the environment types so that we can plot them all on the
    # last plot
    env_type_df = pd.DataFrame.from_items([(s.name, s) for s in series_holder]).T
    # here we need to add the columns to the df.
    # use the sampling frequency values
    env_type_df.columns = sampling_frequencies
    colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                          create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000,
                                             num_cols=len(list(env_type_df)),
                                             time_out_iterations=10000)]
    colour_dict = {type_id: colour for type_id, colour in zip(env_type_df.index.values.tolist(), colour_palette_pas)}
    # now create a colour dict
    type_id_list = env_type_df.index.values.tolist()
    # plot with full sampling, then with only up to 1000
    for type_id in type_id_list:
        # here plot the individual points (one pont for each of the samples of the env_type that have a point for
        # this sampling frequency
        y_list = env_type_df.loc[type_id].values.tolist()
        x_list = list(env_type_df)

        # now draw the line
        axarr[0].plot(x_list, y_list, color=colour_dict[type_id])
    axarr[0].set_xscale('log')

    apples = 'asdf'
# TODO replace the middle chart with the bubble plot with connectivity
def plot_rarefaction_indi_gen_two(type_profile_rel_abund_dict_holder_dict):
    # sampling frequecies
    sampling_frequencies = []
    additions_list = [0, 0.25, 0.5, 0.75]
    orders = range(1, 6)
    for order in orders:
        for addition in additions_list:
            sampling_frequencies.append(int(10 ** (order + addition)))
    if os.path.isfile('type_profile_rarefaction_dict.pickle'):
        output_dict_shared = pickle.load(open('type_profile_rarefaction_dict.pickle', 'rb'))
    else:

        # 1

        # create a dict to put the output in
        worker_manager = Manager()
        output_dict_shared = worker_manager.dict()

        num_proc = 7

        # put in tups that are the items in the type_profile_rel_abund_dict_holder_dict
        input = Queue()
        for tup in type_profile_rel_abund_dict_holder_dict.items():
            input.put(tup)

        for i in range(num_proc):
            input.put('STOP')

        list_of_processes = []
        for N in range(num_proc):
            p = Process(target=rarefaction_curve_worker_indi_type_profile,
                        args=(input, 100, output_dict_shared, sampling_frequencies))
            list_of_processes.append(p)
            p.start()

        for p in list_of_processes:
            p.join()

        pickle.dump(dict(output_dict_shared), open('type_profile_rarefaction_dict.pickle', 'wb'))

    # PLOT RAREFACTION

    fig = plt.figure(figsize=(14, 6))
    gs = plt.GridSpec(1, 6, width_ratios=[3,1.5,1.5,1.5,1.5,3], figure=fig)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])
    ax4 = plt.subplot(gs[4])
    ax5 = plt.subplot(gs[5])

    plt.subplots_adjust(left=0.06, right=0.94, wspace=0.5)
    apples = 'asdf'


    # list that will hold the series that will hold the averaged abundances at each bootstrap for each type profile
    series_holder = []
    # for each typeprofile workout the means of the bottstraps as a series and add these to the list
    for type_prof in output_dict_shared.keys():
        temp_df = pd.DataFrame(output_dict_shared[type_prof])
        averaged_series = temp_df.mean(axis=0)
        averaged_series.name = type_prof
        series_holder.append(averaged_series)
    # here we should have all of the info we need to make a dataframe that can directly be used for plotting
    # This will hold the pairs of x, y lists for each of the environment types so that we can plot them all on the
    # last plot
    line_holder = []
    env_type_df = pd.DataFrame.from_items([(s.name, s) for s in series_holder]).T
    # here we need to add the columns to the df.
    # use the sampling frequency values
    env_type_df.columns = sampling_frequencies
    colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                          create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000,
                                             num_cols=len(list(env_type_df)),
                                             time_out_iterations=10000)]
    colour_dict = {type_id: colour for type_id, colour in zip(env_type_df.index.values.tolist(), colour_palette_pas)}
    # now create a colour dict
    type_id_list = env_type_df.index.values.tolist()

    # holder for the bubble plot info
    dict_100 = {}
    dict_1000 = {}
    dict_10000 = {}
    dict_final = {}

    # plot with full sampling, then with only up to 1000
    for type_id in type_id_list:
        # here plot the individual points (one pont for each of the samples of the env_type that have a point for
        # this sampling frequency
        y_list = env_type_df.loc[type_id].values.tolist()
        x_list = list(env_type_df)

        dict_100[type_id] = y_list[4]
        dict_1000[type_id] = y_list[8]
        dict_10000[type_id] = y_list[12]
        dict_10000[type_id] = y_list[16]
        dict_final[type_id] = y_list[-1]
        # now draw the line
        ax0.plot(x_list, y_list, color='black')
        ax0.text(max(x_list), max(y_list), str(type_id))
    ax0.set_xlim(10, 1000000)
    ax0.set_xscale('symlog')
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)


    # this is where we will want to do the bubble plot things
    # we y coordinates should be got by dividing up the y by num of types + 1 pieces

    #100
    plot_bubble(axes=ax1, predictor_dict=dict_100, final_dict=dict_final, type_id_list=type_id_list, label='100 sequences', colour='#D3D3D3')

    #1000
    plot_bubble(axes=ax2, predictor_dict=dict_1000, final_dict=dict_final, type_id_list=type_id_list, label='1000 sequences', colour='#A9A9A9')

    #10000
    plot_bubble(axes=ax3, predictor_dict=dict_1000, final_dict=dict_final, type_id_list=type_id_list, label='10000 sequences', colour='#696969')

    # 100000
    plot_bubble(axes=ax4, predictor_dict=dict_10000, final_dict=dict_final, type_id_list=type_id_list,
                label='100000 sequences', colour='black')

    # now make a scatter plot with a regresion line for the 100, 1 000, 10 000 sampling
    # showing the value vs the known end
    colour_dict_sampling = {100: '#D3D3D3', 1000: "#A9A9A9", 10000: "#696969", 100000: 'black'}
    for sampling in [100, 1000, 10000, 100000]:
        x_list = env_type_df[sampling].values.tolist()
        y_list = env_type_df[562341].values.tolist()
        # now draw the line
        ax5.scatter(x_list, y_list, color=colour_dict_sampling[sampling], marker='.')
        # https://stackoverflow.com/questions/22239691/code-for-line-of-best-fit-of-a-scatter-plot-in-python
        ax5.plot(np.unique(x_list), np.poly1d(np.polyfit(x_list, y_list, 1))(np.unique(x_list)), color=colour_dict_sampling[sampling])

    ax5.set_xlim(0, 1200)
    ax5.set_ylim(0, 1200)
    # now plot the 1:1 line

    # ax5.plot(np.unique(x_list), np.unique(x_list), linestyle='--', color='black')

    # plot the 1 for 1 line on the third subplot so that we can see what a good predictor it is
    plt.savefig('type_profile_individual_rarefaction.svg')
    plt.savefig('type_profile_individual_rarefaction.png')


def plot_bubble(axes, predictor_dict, final_dict, type_id_list, label, colour):
    axes.set_xlim(-1.2, 2.2)
    axes.set_ylim(0, 1)
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off',
                    labelbottom='off', labeltop='off', labelleft='off', labelright='off', )
    axes.arrow(-1, 0, 0, 0.99, head_width=0.2, head_length=0.05, length_includes_head=True, color='gray',
              edgecolor='gray')
    axes.arrow(2, 0, 0, 0.99, head_width=0.2, head_length=0.05, length_includes_head=True, color='black',
              edgecolor='black')
    # axes.text(-1, 1.02, label)
    axes.set_title(label)
    num_types = len(type_id_list)
    # patch list to hold the circles.
    patches_list = []
    y_cords_list = [0 + (i * (1 / (num_types + 1))) for i in range(num_types + 1)][1:]
    x_coords = [-0.5, 1.5]
    # now we simply put the points into order
    # first plot the predictor order
    # we should have a dict of the coordinate for each type for making the line
    predictor_coord_dict = {}
    end_abund_coord_dict = {}
    sorted_predictors = sorted(predictor_dict.items(), key=lambda x: x[1])
    sorted_end_abund = sorted(final_dict.items(), key=lambda x: x[1])
    # I didn't use the circle path in the end becuase getting the true circle was a pain with unequal axes
    for i, tup in enumerate(sorted_predictors):
        x_y = (x_coords[0], y_cords_list[i])
        axes.scatter(x_coords[0], y_cords_list[i], marker='o', s=100, color=colour, zorder=2)
        # patches_list.append(Circle(xy=x_y, radius=0.5))
        predictor_coord_dict[tup[0]] = x_y
    for i, tup in enumerate(sorted_end_abund):
        x_y = (x_coords[1], y_cords_list[i])
        axes.scatter(x_coords[1], y_cords_list[i], marker='o', s=100, facecolors='none', edgecolors='black', zorder=2)
        # patches_list.append(Circle(xy=x_y, radius=0.2))
        end_abund_coord_dict[tup[0]] = x_y
    # patches_collection = PatchCollection(patches_list)
    # ax1.add_collection(patches_collection)
    # now draw the lines between the two points
    for type_id in predictor_coord_dict.keys():
        x_s = [predictor_coord_dict[type_id][0], end_abund_coord_dict[type_id][0]]
        y_s = [predictor_coord_dict[type_id][1], end_abund_coord_dict[type_id][1]]
        axes.plot(x_s, y_s, c='black', zorder=1)


def rarefaction_curve_worker_indi_type_profile(input_queue, num_bootstraps, result_dict, sampling_frequencies):

    for name, working_dict in iter(input_queue.get, 'STOP'):
        fixed_dict_tup = list(working_dict.items())
        sequences = [tup[0] for tup in fixed_dict_tup]
        probabilities = [tup[1] for tup in fixed_dict_tup]
        # hardcode check that the probabilities == 1
        if sum(probabilities) != 1:
            if sum(probabilities) > 1:
                probabilities[0] -= sum(probabilities) - 1
            else:
                probabilities[0] += 1 - sum(probabilities)

        sys.stdout.write('\n\nType: {}. Probs sum to {}'.format(name, sum(probabilities)))
        # for each sample, perform the bootstrapping

        # this will hold the lists which will hold the results of a single bool
        # so this will hold the results of all of the bools
        sample_boot_result_holder = []
        for it_num in range(num_bootstraps):
            sys.stdout.write('\nbootstrap: {}\n'.format(it_num))
            # this is the set that we will populate to count number of unique seqs
            sample_result_holder = []

            picking_list = np.random.choice(a=sequences, size=max(sampling_frequencies), p=probabilities)

            for i in sampling_frequencies:
                sys.stdout.write('\rsampling at: {}'.format(i))
                sample_result_holder.append(len(set(picking_list[:i])))

            sample_boot_result_holder.append(sample_result_holder)

        # here we have conducted the bootstrapping for one of the samples
        result_dict[name] = sample_boot_result_holder

def evenness_worker_indi_type_profile(input_queue, num_bootstraps, result_dict, sampling_frequencies):

    for name, working_dict in iter(input_queue.get, 'STOP'):
        fixed_dict_tup = list(working_dict.items())
        sequences = [tup[0] for tup in fixed_dict_tup]
        probabilities = [tup[1] for tup in fixed_dict_tup]
        # hardcode check that the probabilities == 1
        if sum(probabilities) != 1:
            if sum(probabilities) > 1:
                probabilities[0] -= sum(probabilities) - 1
            else:
                probabilities[0] += 1 - sum(probabilities)

        sys.stdout.write('\n\nType: {}. Probs sum to {}'.format(name, sum(probabilities)))
        # for each sample, perform the bootstrapping

        # this will hold the lists which will hold the results of a single bool
        # so this will hold the results of all of the bools
        sample_boot_result_holder = []
        for it_num in range(num_bootstraps):
            sys.stdout.write('\nbootstrap: {}\n'.format(it_num))
            # this is the set that we will populate to count number of unique seqs
            sample_result_holder = []

            picking_list = np.random.choice(a=sequences, size=max(sampling_frequencies), p=probabilities)

            for i in sampling_frequencies:
                sys.stdout.write('\rsampling at: {}'.format(i))
                sample_result_holder.append(calculate_shannons_equitability(picking_list[:i]))

            sample_boot_result_holder.append(sample_result_holder)

        # here we have conducted the bootstrapping for one of the samples
        result_dict[name] = sample_boot_result_holder



def create_colour_list(sq_dist_cutoff=None, mix_col=None, num_cols=50, time_out_iterations=10000, avoid_black_and_white=True):
    new_colours = []
    min_dist = []
    attempt = 0
    while len(new_colours) < num_cols:
        attempt += 1
        # Check to see if we have run out of iteration attempts to find a colour that fits into the colour space
        if attempt > time_out_iterations:
            sys.exit('Colour generation timed out. We have tried {} iterations of colour generation '
                     'and have not been able to find a colour that fits into your defined colour space.\n'
                     'Please lower the number of colours you are trying to find, '
                     'the minimum distance between them, or both.'.format(attempt))
        if mix_col:
            r = int((random.randint(0, 255) + mix_col[0]) /2)
            g = int((random.randint(0, 255) + mix_col[1]) /2)
            b = int((random.randint(0, 255) + mix_col[2]) /2)
        else:
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        # now check to see whether the new colour is within a given distance
        # if the avoids are true also
        good_dist = True
        if sq_dist_cutoff:
            dist_list = []
            for i in range(len(new_colours)):
                distance = (new_colours[i][0] - r)**2 + (new_colours[i][1] - g)**2 + (new_colours[i][2] - b)**2
                dist_list.append(distance)
                if distance < sq_dist_cutoff:
                    good_dist = False
                    break
            # now check against black and white
            d_to_black = (r - 0)**2 + (g - 0)**2 + (b - 0)**2
            d_to_white = (r - 255)**2 + (g - 255)**2 + (b - 255)**2
            if avoid_black_and_white:
                if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                    good_dist = False
            if dist_list:
                min_dist.append(min(dist_list))
        if good_dist:
            new_colours.append((r,g,b))
            attempt = 0

    return new_colours

def calculate_shannons_equitability(list_of_items):
    ''' This method will take in a list of items. It will calculate the shannon's eveness for this list.
    To do this it will first create a dictionary that is a counter. It will then use this to do
    the caluclation.'''

    temp_counting_dict = Counter(list_of_items)

    # convert the counter to a plain old dict
    temp_counting_dict = dict(temp_counting_dict)

    # get total number of sequences
    total_number_of_seqs = sum(list(temp_counting_dict.values()))

    # get number of distinct seqs
    distinct_seqs = len(temp_counting_dict.keys())

    # calculate shannon
    shannon = 0
    for val in temp_counting_dict.values():
        proportion = val/total_number_of_seqs
        natural_log_of_proportion = math.log(proportion)
        shannon += proportion * natural_log_of_proportion

    shannon = shannon * -1

    shannons_evenness = shannon / math.log(distinct_seqs)

    return shannons_evenness



def generate_median_joining_network_from_dict(seq_abund_dict, type_id, subsample_level=1000):
    '''

    I am going to modify this so that it also produces a MED version for plotting up.
    The aim of this method will be to generate a network for a set of sequences. To start with this will be
    passed in as a dict with sequence as key and abundance as value. We will work with these abundances as relative
    so that the largest is 1 and the smallest will be scaled relative.
    1 - Blast the sequences against a clade db to get clade. Infer the major clade of the type profile
    then dicard sequeneces that are not from the same clade.
    2 - produce an alignment of the sequences
    3 - put this alignment into splits trees somehow
    4 - transfer this into network x and make a network from it'''

    # it may be best if we pass in a redundant fasta to splits tree as then the actual sizes will be correct
    # then perhaps we can even draw this ourselves using circles and lines and get the sizes directly from the
    # splitstree output file

    # create a fasta file from the dictionary
    # we will work with a set size of 100000 sequences
    type_wkd = '{}/{}'.format(os.getcwd(), type_id)
    os.makedirs(type_wkd, exist_ok=True)

    if os.path.isfile(
            '{}/{}_splits_tree_out_file_subsampled_no_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level)) and os.path.isfile(
            '{}/{}_splits_tree_out_file_subsampled_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level)):

        splits_tree_out_file_no_med = pickle.load(open('{}/{}_splits_tree_out_file_subsampled_no_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level), 'rb'))
        count_id_to_abund_dict = pickle.load(open('{}/{}_count_id_to_abund_dict_subsampled_no_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level), 'rb'))
        splits_tree_out_file_med = pickle.load(open('{}/{}_splits_tree_out_file_subsampled_med_{}_dynamicM.pickle'.format(type_wkd,type_id, subsample_level), 'rb'))
        med_abundance_dict = pickle.load(open('{}/{}_count_id_to_abund_dict_subsampled_med_{}_dynamicM.pickle'.format(type_wkd,type_id, subsample_level), 'rb'))

    else:

        fasta_to_blast = []
        count = 0
        for sequence, abundance in seq_abund_dict.items():
            # here check to see that the sequence will be found still when subsampling to 1000
            if int(abundance*subsample_level) > 0:
                fasta_to_blast.extend(['>seq{}'.format(count), '{}'.format(sequence)])
            count += 1

        # create a dictionary of name to seq
        count_id_to_seq_dict = {fasta_to_blast[i][1:] : fasta_to_blast[i+1] for i in range(0, len(fasta_to_blast), 2)}

        writeListToDestination('{}/seqs_to_write.fasta'.format(os.getcwd()), fasta_to_blast)

        # now perform the blast and get a sequence to clade dictionary as a result
        # now do a blast on these seqs
        ncbircFile = []
        if os.path.isdir('/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/symbiodiniumDB'):
            db_path = '/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/symbiodiniumDB'
        else:
            # if we are on the remote server
            db_path = '/home/humebc/phylogeneticSoftware/SymPortal_framework/symbiodiniumDB'
        ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])
        # write the .ncbirc file that gives the location of the db
        writeListToDestination('{}/.ncbirc'.format(os.getcwd()), ncbircFile)
        blastOutputPath = '{}/blast.out'.format(os.getcwd())
        outputFmt = "6 qseqid sseqid staxids evalue"
        inputPath = '{}/seqs_to_write.fasta'.format(os.getcwd())

        # Run local blast
        # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
        completedProcess = subprocess.run(
            ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db',
             'symClade.fa', '-max_target_seqs', '1', '-num_threads', '7'])

        # Read in blast output
        blast_output_file = readDefinedFileToList(blastOutputPath)

        seq_clade_dict = {}
        for result_line in blast_output_file:
            seq_clade_dict[result_line.split('\t')[0]] = result_line.split('\t')[1][-1]

        # get the most common clade
        clade_counter = defaultdict(int)
        for value in seq_clade_dict.values():
            clade_counter[value] += 1

        maj_clade = sorted(clade_counter.items(), key=lambda x: x[1], reverse=True)[0][0]

        # now discard sequences from the dictionary that do not match the maj_clade
        non_clade_seq_counter = 0
        for seq, clade in seq_clade_dict.items():
            if clade != maj_clade:
                seq_to_del = count_id_to_seq_dict[seq]
                del seq_abund_dict[seq_to_del]
                non_clade_seq_counter +=1
        print('{} sequences deleted'.format(non_clade_seq_counter))

        # now recreate the fasta from the dict with only single clade seqs
        # create a fasta file from the dictionary
        fasta_to_align = []
        count = 0
        for sequence, abundance in seq_abund_dict.items():
            if int(abundance * subsample_level) > 0:
                fasta_to_align.extend(['>seq{}'.format(count), '{}'.format(sequence)])
                count += 1


        # we need to have a count_id_to_abundance_dict for drawing the networks
        count_id_to_abund_dict = {fasta_to_align[i][1:]: seq_abund_dict[fasta_to_align[i + 1]] for i in range(0, len(fasta_to_align), 2)}

        # here we have a cleaned dict that only contains sequences from one clade
        # now align using mafft
        # Write out the new fasta
        infile_path = '{}/temp_type_fasta_for_net_blast.fasta'.format(os.getcwd())
        writeListToDestination(infile_path, fasta_to_align)
        # now perform the alignment with MAFFT
        mafft = local["mafft-linsi"]

        out_file = infile_path.replace('.fasta', '_aligned.fasta')

        # now run mafft including the redirect
        (mafft[ '--thread', -1, infile_path] > out_file)()

        aligned_fasta_interleaved = readDefinedFileToList(out_file)
        aligned_fasta = convert_interleaved_to_sequencial_fasta_two(aligned_fasta_interleaved)
        aligned_fasta_cropped = crop_fasta(aligned_fasta)


        # cropped_count_id_to_abundance_dict = {}
        # for i in range(0, len(aligned_fasta_cropped), 2):
        #     count_id = aligned_fasta_cropped[i][1:]
        #     original_seq = count_id_to_seq_dict[count_id]
        #     abundance = seq_abund_dict[original_seq]
        #     cropped_count_id_to_abundance_dict[count_id] = abundance

        # here we now have the fasta cropped and a new copped_seq_abundance_dict
        # we also have a dictionary that links the count_id to the abundance which is what we can use
        # when creating the network

        # here is where we should implement MED
        # it would be good to have a method that took in the aligned_fasta, the dict and then returned a new fasta
        # and dict
        med_aligned_fasta, med_abundance_dict = perform_MED_on_fasta(aligned_fasta_cropped=aligned_fasta_cropped,
                                                                     count_id_to_abund=count_id_to_abund_dict,
                                                                     subsample=subsample_level, type_id=type_id)

        # we have to make the new nexus format by hand as the biopython version was putting out old versions.
        new_nexus_no_med = splits_tree_nexus_from_fasta(aligned_fasta_cropped)
        new_nexus_med = splits_tree_nexus_from_fasta(med_aligned_fasta)

        # now the nexus is complete and ready for writing out.
        # I think we whould start working in directories here
        network_wkd = '{}/{}/{}_networks'.format(os.getcwd(), type_id, type_id)
        os.makedirs(network_wkd, exist_ok=True)

        splits_nexus_path_no_med = '{}/{}_splitstree_in_no_med.nex'.format(network_wkd, type_id)
        with open(splits_nexus_path_no_med, 'w') as f:
            for line in new_nexus_no_med:
                f.write('{}\n'.format(line))

        splits_nexus_path_med = '{}/{}_splitstree_in_med.nex'.format(network_wkd, type_id)
        with open(splits_nexus_path_med, 'w') as f:
            for line in new_nexus_med:
                f.write('{}\n'.format(line))


        # now create the control file that we can use for execution for the no med
        splits_out_path_no_med = '{}/{}_splitstree_out_no_med.nex'.format(network_wkd, type_id)
        ctrl_out_path_no_med = '{}/{}_splitstree_ctrl_no_med'.format(network_wkd, type_id)
        # this creates a control file that can be fed to splits tree on the command line and write it out
        # it then tuns splitstrees with the cntrl file before returning the output file
        splits_tree_out_file_no_med = run_splits_trees(ctrl_out_path=ctrl_out_path_no_med, splits_nexus_path=splits_nexus_path_no_med, splits_out_path=splits_out_path_no_med)

        # now create the control file that we can use for execution for the no med
        splits_out_path_med = '{}/{}_splitstree_out_med.nex'.format(network_wkd, type_id)
        ctrl_out_path_med = '{}/{}_splitstree_ctrl_med'.format(network_wkd, type_id)
        # this creates a control file that can be fed to splits tree on the command line and write it out
        splits_tree_out_file_med = run_splits_trees(ctrl_out_path=ctrl_out_path_med,
                                                    splits_nexus_path=splits_nexus_path_med,
                                                    splits_out_path=splits_out_path_med)


        #TODO perhaps we can just pass out the count_id_to abundance dict rather than both
        # this should save us a step later on. Do the same for the MED version of this.
        pickle.dump(splits_tree_out_file_no_med, open(
            '{}/{}_splits_tree_out_file_subsampled_no_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level), 'wb'))
        pickle.dump(count_id_to_abund_dict, open(
            '{}/{}_count_id_to_abund_dict_subsampled_no_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level),
            'wb'))
        pickle.dump(splits_tree_out_file_med, open(
            '{}/{}_splits_tree_out_file_subsampled_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level), 'wb'))
        pickle.dump(med_abundance_dict, open(
            '{}/{}_count_id_to_abund_dict_subsampled_med_{}_dynamicM.pickle'.format(type_wkd, type_id, subsample_level), 'wb'))



    return [(splits_tree_out_file_no_med, count_id_to_abund_dict), (splits_tree_out_file_med, med_abundance_dict)]

def generate_median_joining_network_from_dict_no_med(seq_abund_dict, type_id, subsample_level=1000):
    '''

    I am going to modify this so that it also produces a MED version for plotting up.
    The aim of this method will be to generate a network for a set of sequences. To start with this will be
    passed in as a dict with sequence as key and abundance as value. We will work with these abundances as relative
    so that the largest is 1 and the smallest will be scaled relative.
    1 - Blast the sequences against a clade db to get clade. Infer the major clade of the type profile
    then dicard sequeneces that are not from the same clade.
    2 - produce an alignment of the sequences
    3 - put this alignment into splits trees somehow
    4 - transfer this into network x and make a network from it'''

    # it may be best if we pass in a redundant fasta to splits tree as then the actual sizes will be correct
    # then perhaps we can even draw this ourselves using circles and lines and get the sizes directly from the
    # splitstree output file

    # create a fasta file from the dictionary
    # we will work with a set size of 100000 sequences
    type_wkd = '{}/{}'.format(os.getcwd(), type_id)
    os.makedirs(type_wkd, exist_ok=True)

    if os.path.isfile(
            '{}/{}_splits_tree_out_file_subsampled_no_med_{}.pickle'.format(type_wkd, type_id, subsample_level)):

        splits_tree_out_file_no_med = pickle.load(open('{}/{}_splits_tree_out_file_subsampled_no_med_{}.pickle'.format(type_wkd, type_id, subsample_level), 'rb'))
        count_id_to_abund_dict = pickle.load(open('{}/{}_count_id_to_abund_dict_subsampled_no_med_{}.pickle'.format(type_wkd, type_id, subsample_level), 'rb'))

    else:
        fasta_to_blast = []
        count = 0
        for sequence, abundance in seq_abund_dict.items():
            # here check to see that the sequence will be found still when subsampling to 1000
            if int(abundance*subsample_level) > 0:
                fasta_to_blast.extend(['>seq{}'.format(count), '{}'.format(sequence)])
            count += 1

        # create a dictionary of name to seq
        count_id_to_seq_dict = {fasta_to_blast[i][1:] : fasta_to_blast[i+1] for i in range(0, len(fasta_to_blast), 2)}

        writeListToDestination('{}/seqs_to_write.fasta'.format(os.getcwd()), fasta_to_blast)

        # now perform the blast and get a sequence to clade dictionary as a result
        # now do a blast on these seqs
        ncbircFile = []
        if os.path.isdir('/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/symbiodiniumDB'):
            db_path = '/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/symbiodiniumDB'
        else:
            # if we are on the remote server
            db_path = '/home/humebc/phylogeneticSoftware/SymPortal_framework/symbiodiniumDB'
        ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])
        # write the .ncbirc file that gives the location of the db
        writeListToDestination('{}/.ncbirc'.format(os.getcwd()), ncbircFile)
        blastOutputPath = '{}/blast.out'.format(os.getcwd())
        outputFmt = "6 qseqid sseqid staxids evalue"
        inputPath = '{}/seqs_to_write.fasta'.format(os.getcwd())

        # Run local blast
        # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
        completedProcess = subprocess.run(
            ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db',
             'symClade.fa', '-max_target_seqs', '1', '-num_threads', '7'])

        # Read in blast output
        blast_output_file = readDefinedFileToList(blastOutputPath)

        seq_clade_dict = {}
        for result_line in blast_output_file:
            seq_clade_dict[result_line.split('\t')[0]] = result_line.split('\t')[1][-1]

        # get the most common clade
        clade_counter = defaultdict(int)
        for value in seq_clade_dict.values():
            clade_counter[value] += 1

        maj_clade = sorted(clade_counter.items(), key=lambda x: x[1], reverse=True)[0][0]

        # now discard sequences from the dictionary that do not match the maj_clade
        non_clade_seq_counter = 0
        for seq, clade in seq_clade_dict.items():
            if clade != maj_clade:
                seq_to_del = count_id_to_seq_dict[seq]
                del seq_abund_dict[seq_to_del]
                non_clade_seq_counter +=1
        print('{} sequences deleted'.format(non_clade_seq_counter))

        # now recreate the fasta from the dict with only single clade seqs
        # create a fasta file from the dictionary
        fasta_to_align = []
        count = 0
        for sequence, abundance in seq_abund_dict.items():
            if int(abundance * subsample_level) > 0:
                fasta_to_align.extend(['>seq{}'.format(count), '{}'.format(sequence)])
                count += 1


        # we need to have a count_id_to_abundance_dict for drawing the networks
        count_id_to_abund_dict = {fasta_to_align[i][1:]: seq_abund_dict[fasta_to_align[i + 1]] for i in range(0, len(fasta_to_align), 2)}

        # here we have a cleaned dict that only contains sequences from one clade
        # now align using mafft
        # Write out the new fasta
        infile_path = '{}/temp_type_fasta_for_net_blast.fasta'.format(os.getcwd())
        writeListToDestination(infile_path, fasta_to_align)
        # now perform the alignment with MAFFT
        mafft = local["mafft-linsi"]

        out_file = infile_path.replace('.fasta', '_aligned.fasta')

        # now run mafft including the redirect
        (mafft[ '--thread', -1, infile_path] > out_file)()

        aligned_fasta_interleaved = readDefinedFileToList(out_file)
        aligned_fasta = convert_interleaved_to_sequencial_fasta_two(aligned_fasta_interleaved)
        aligned_fasta_cropped = crop_fasta(aligned_fasta)




        # here we now have the fasta cropped and a new copped_seq_abundance_dict
        # we also have a dictionary that links the count_id to the abundance which is what we can use
        # when creating the network

        # we have to make the new nexus format by hand as the biopython version was putting out old versions.
        new_nexus_no_med = splits_tree_nexus_from_fasta(aligned_fasta_cropped)

        # now the nexus is complete and ready for writing out.
        # I think we whould start working in directories here
        network_wkd = '{}/{}/{}_networks'.format(os.getcwd(), type_id, type_id)
        os.makedirs(network_wkd, exist_ok=True)

        splits_nexus_path_no_med = '{}/{}_splitstree_in_no_med.nex'.format(network_wkd, type_id)
        with open(splits_nexus_path_no_med, 'w') as f:
            for line in new_nexus_no_med:
                f.write('{}\n'.format(line))




        # now create the control file that we can use for execution for the no med
        splits_out_path_no_med = '{}/{}_splitstree_out_no_med.nex'.format(network_wkd, type_id)
        ctrl_out_path_no_med = '{}/{}_splitstree_ctrl_no_med'.format(network_wkd, type_id)
        # this creates a control file that can be fed to splits tree on the command line and write it out
        # it then tuns splitstrees with the cntrl file before returning the output file
        splits_tree_out_file_no_med = run_splits_trees(ctrl_out_path=ctrl_out_path_no_med, splits_nexus_path=splits_nexus_path_no_med, splits_out_path=splits_out_path_no_med)



        #TODO perhaps we can just pass out the count_id_to abundance dict rather than both
        # this should save us a step later on. Do the same for the MED version of this.
        pickle.dump(splits_tree_out_file_no_med, open(
            '{}/{}_splits_tree_out_file_subsampled_no_med_{}.pickle'.format(type_wkd, type_id, subsample_level), 'wb'))
        pickle.dump(count_id_to_abund_dict, open(
            '{}/{}_count_id_to_abund_dict_subsampled_no_med_{}.pickle'.format(type_wkd, type_id, subsample_level),
            'wb'))



    return splits_tree_out_file_no_med, count_id_to_abund_dict


def run_splits_trees(ctrl_out_path, splits_nexus_path, splits_out_path):
    ctrl_file = []
    ctrl_file.append('BEGIN SplitsTree;')
    ctrl_file.append('EXECUTE FILE={};'.format(splits_nexus_path))
    ctrl_file.append('SAVE FILE={} REPLACE=yes;'.format(splits_out_path))
    ctrl_file.append('QUIT;')
    ctrl_file.append('end;')
    # now write out the control file
    with open(ctrl_out_path, 'w') as f:
        for line in ctrl_file:
            f.write('{}\n'.format(line))

    # now run splitstree
    completedProcess = subprocess.run(
        ['SplitsTree', '-g', '-c', ctrl_out_path])

    # now we can read in the output file
    # and then we can start making the network finally!
    with open(splits_out_path, 'r') as f:
        splits_tree_out_file = [line.rstrip() for line in f]

    return splits_tree_out_file

def crop_fasta(aligned_fasta):
    # convert each of the sequences in the fasta into a series with the series name as the sequence name from the fasta
    temp_series_list = []
    for i in range(0, len(aligned_fasta), 2):
        temp_series_list.append(pd.Series(list(aligned_fasta[i+1]), name=aligned_fasta[i][1:]))

    # now create the df from the list of series
    # https://github.com/pandas-dev/pandas/issues/1494
    aligned_fasta_as_df = pd.DataFrame.from_items([(s.name, s) for s in temp_series_list]).T
    # aligned_fasta_as_df = pd.DataFrame(temp_series_list)

    # now do the cropping
    aligned_fasta_as_df_cropped = crop_fasta_df(aligned_fasta_as_df)

    # now we need to convert this back to a fasta
    output_fasta = []
    for sequence in aligned_fasta_as_df_cropped.index.values.tolist():
        output_fasta.extend(['>{}'.format(aligned_fasta_as_df_cropped.loc[sequence].name), ''.join(aligned_fasta_as_df_cropped.loc[sequence].values.tolist())])

    return output_fasta

def crop_fasta_df(aligned_fasta_as_pandas_df_to_crop):
    columns_to_drop = []
    for i in list(aligned_fasta_as_pandas_df_to_crop):
        # if there is a gap in the column at the beginning
        if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
            columns_to_drop.append(i)
        else:
            break
    for i in reversed(list(aligned_fasta_as_pandas_df_to_crop)):
        # if there is a gap in the column at the end
        if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
            columns_to_drop.append(i)
        else:
            break

    # get a list that is the columns indices that we want to keep
    col_to_keep = [col_index for col_index in list(aligned_fasta_as_pandas_df_to_crop) if col_index not in columns_to_drop]
    # drop the gap columns
    return aligned_fasta_as_pandas_df_to_crop[col_to_keep]

def perform_MED_on_fasta(aligned_fasta_cropped, count_id_to_abund, subsample, type_id):
    # so the basic idea here is to create a redundant fasta that we can write out and then feed into the MED
    # decomposition. Once we have done the MED we will have to work through the output to get which sequences that are the MED nodes
    # we will then need to get the abundances and use the combination of these two to make a dictionary
    # This dictionary will be returned. The fasta that we put into the MED should contain unaligned sequences.
    # these sequences will need to be padded before we run the MED.

    sys.stdout.write('\n{}: padding alignment\n'.format(type_id))

    # convert to an unaligned fasta
    unaligned_fasta = remove_gaps_from_aligned_fasta(aligned_fasta_cropped)

    MEDOutDir = '{}/{}/{}_med_out'.format(os.getcwd(), type_id, type_id)
    os.makedirs(MEDOutDir, exist_ok=True)

    # write out the fasta so that it can be padded and then run through MED
    path_to_unalined_fasta = '{}/{}_cropped_fasta_for_padding.fasta'.format(MEDOutDir, type_id)
    with open(path_to_unalined_fasta, 'w') as f:
        for line in unaligned_fasta:
            f.write('{}\n'.format(line))

    completedProcess = subprocess.run([r'o-pad-with-gaps', r'{}'.format(path_to_unalined_fasta)])

    # now we need to infer the path of the padded files and then pass these into the med
    # TODO do this by looking
    padded_file_path = path_to_unalined_fasta + '-PADDED-WITH-GAPS'

    # now we need to read this back in so that we can make the fasta redundant so that MED can do its thing
    with open(padded_file_path, 'r') as f:
        padded_fasta = [line.rstrip() for line in f]

    # multioply out the relative abundances by the subsampling level
    reduntant_padded_fasta = []
    count = 0
    for i in range(0, len(padded_fasta), 2):
        relabundance = count_id_to_abund[padded_fasta[i][1:]]
        sub_sampled_abundance = int(subsample * relabundance)
        if sub_sampled_abundance > 0:
            for j in range(sub_sampled_abundance):
                reduntant_padded_fasta.append('>seq{}'.format(count))
                reduntant_padded_fasta.append(padded_fasta[i+1])
                count += 1

    # now we can write out the padded fasta and feed this into the med
    path_to_cropped_padded_redundant_fasta = '{}/{}_cropped_fasta_padded_redundant.fasta'.format(MEDOutDir, type_id)
    with open(path_to_cropped_padded_redundant_fasta, 'w') as f:
        for line in reduntant_padded_fasta:
            f.write('{}\n'.format(line))


    sys.stdout.write('{}: running MED\n'.format(type_id))

    #Here we need to make sure that the M value is set dynamically. We want to maintain the 1000 to 4 ratio
    # i.e. we want it set to 4 when we are subsampling to 1000. So, we need to set it to 0.04% or 0.004
    M_value = M_value = max(4, int(0.004 * subsample))
    completedProcess = subprocess.run(
        [r'decompose', '-M', str(M_value), '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html', '--skip-check-input', '-o',
         MEDOutDir, path_to_cropped_padded_redundant_fasta])

    # we can work with the NODE-REPRESENTATIVES.fasta file to read in the representative sequence. This
    # also has abundance information with it.
    # We can have the opportunity to really play around with the effect of the M value here. For the time
    # being we will continue to work with -M 1 but we should certainly take a look at how the deafault -M 4 performs
    # one thing I would like to check is to see that the umber of samples in the NODE-REp file actually
    # add up to what we are expecting.
    with open('{}/NODE-REPRESENTATIVES.fasta'.format(MEDOutDir), 'r') as f:
        node_rep_file = [line.rstrip() for line in f]

    med_fasta_unaligned = []
    count = 0

    # we want to return the relative abundance rather than the med absolute abundances
    # so divide the MED absolute abundances for each node by the total
    med_count_id_to_abundance = {}
    total_med_node_abundance = sum([int(node_rep_file[i].split(':')[-1]) for i in range(0, len(node_rep_file), 2)])
    for i in range(0, len(node_rep_file), 2):
        med_fasta_unaligned.append('>seq{}'.format(count))
        med_fasta_unaligned.append(node_rep_file[i+1].replace('-', ''))
        med_count_id_to_abundance['seq{}'.format(count)] = float(int(node_rep_file[i].split(':')[-1]) / total_med_node_abundance)
        count += 1

    # write out the unaligned med fasta
    unaligned_post_med_fasta_path = '{}/{}_unaligned_post_med.fasta'.format(MEDOutDir, type_id)
    with open(unaligned_post_med_fasta_path, 'w') as f:
        for line in med_fasta_unaligned:
            f.write('{}\n'.format(line))


    # now we need to re-align the sequences before we can make a network from them


    # now perform the alignment with MAFFT
    mafft = local["mafft-linsi"]

    out_file = '{}/{}_post_med_aligned.fasta'.format(MEDOutDir, type_id)

    # now run mafft including the redirect
    (mafft['--thread', -1, unaligned_post_med_fasta_path] > out_file)()

    aligned_post_med_fasta_interleaved = readDefinedFileToList(out_file)
    aligned_post_med_fasta = convert_interleaved_to_sequencial_fasta_two(aligned_post_med_fasta_interleaved)
    aligned_post_med_fasta_cropped = crop_fasta(aligned_post_med_fasta)

    # at this point we have a fasta alignment that has been through med, and has been cropped
    # we also have the med_count_id_to_abundance which will work to keep track of the abundances for the network
    # calculation. We should return these both.

    sys.stdout.write('{}: MED complete\n'.format(type_id))

    # read in the med files.


    return aligned_post_med_fasta_cropped, med_count_id_to_abundance

def remove_gaps_from_aligned_fasta(aligned_fasta_list):
    unaligned_fasta = []
    for i in range(0, len(aligned_fasta_list), 2):
        unaligned_fasta.append(aligned_fasta_list[i])
        unaligned_fasta.append((aligned_fasta_list[i+1].replace('-', '')))

    return unaligned_fasta

def splits_tree_nexus_from_fasta(aligned_fasta):
    new_nexus = []
    new_nexus.append('#NEXUS')
    new_nexus.append('BEGIN taxa;')
    new_nexus.append('\tDIMENSIONS ntax={};'.format(int(len(aligned_fasta) / 2)))
    new_nexus.append('TAXLABELS')
    count = 1
    for i in range(0, len(aligned_fasta), 2):
        new_nexus.append('[{}]\t{}'.format(count, aligned_fasta[i][1:]))
        count += 1
    new_nexus.append(';')
    new_nexus.append('END;')
    new_nexus.append('BEGIN characters;')
    new_nexus.append('\tDIMENSIONS nchar={};'.format(len(aligned_fasta[1])))
    new_nexus.append('\tFORMAT')
    new_nexus.append('\t\tdatatype=DNA')
    new_nexus.append('\t\tmissing=?')
    new_nexus.append('\t\tgap=-')
    new_nexus.append('\t\tsymbols="A C G T"')
    new_nexus.append('\t\tlabels=left')
    new_nexus.append('\t\tinterleave=no')
    new_nexus.append('\t;')
    new_nexus.append('MATRIX')
    # now put in the sequences in the style of the mammals.nex example from SplitsTree
    for i in range(0, len(aligned_fasta), 2):
        new_nexus.append('{}\t{}'.format(aligned_fasta[i][1:], aligned_fasta[i + 1].upper()))
    new_nexus.append(';')
    new_nexus.append('END;')
    # finally write in the st_assumption block that will tell SplitsTree to calculate the network
    new_nexus.append('BEGIN st_assumptions;')
    new_nexus.append(
        'CHARTRANSFORM=MedianJoining Epsilon=0 SpringEmbedderIterations=1000 LabelEdges=false ShowHaplotypes=false SubdivideEdges=false ScaleNodesByTaxa=true;')
    new_nexus.append('end;')
    return new_nexus

def draw_network(splits_tree_out_file, count_id_to_abund_dict, ax, colour_for_network, type_id):
    # networkx graph object
    g = nx.Graph()



    # we can work with a minimised version of the splits out file that is only the network block
    for i in range(len(splits_tree_out_file)):
        if splits_tree_out_file[i] == 'BEGIN Network;':
            network_block = splits_tree_out_file[i:]
            break

    # here we can now work with the network_block rather than the larger splits_tree_out_file
    # get a list of the nodes
    # these can be got from the splits_tree_outfile
    # we should use the numbered system that the splits tree is using. We can get a list of the different sequences
    # that each of the taxa represent form the translate part of the network block.
    # This way we can use this to run the sequences against SP to get the names of the sequencs
    # we can also make the 'translation dictionary' at the same time
    vertice_id_to_seq_id_dict = defaultdict(list)
    for i in range(len(network_block)):
        if network_block[i] == 'TRANSLATE':
            # for each line in the translate section
            for j in range(i + 1, len(network_block)):
                if network_block[j] == ';':
                    # then we have reached the end of the translation section
                    break
                items = network_block[j].replace('\'', '').replace(',', '').split(' ')
                vertice_id_to_seq_id_dict[items[0]].extend(items[1:])

    vertices = list(vertice_id_to_seq_id_dict.keys())

    # here we will use the dictionary that was created in the previous method (passed into this one)
    # to link the count id to the actual sequence, which can then be used to look up the abundance
    # This way we can have a vertice_id_to_rel_abund dict that we can use to put a size to the nodes

    vertice_id_to_rel_abund_dict = {}
    for vert_id in vertices:
        count_id_list = vertice_id_to_seq_id_dict[vert_id]
        if len(count_id_list) == 1:
            rel_abund = count_id_to_abund_dict[count_id_list[0]]
            vertice_id_to_rel_abund_dict[vert_id] = rel_abund
        elif len(count_id_list) > 1:
            # then we need sum the abundances for each of them
            tot = 0
            for count_id in count_id_list:
                rel_abund = count_id_to_abund_dict[count_id]
                tot += rel_abund
            vertice_id_to_rel_abund_dict[vert_id] = tot


    # the sizes are a little tricky to work out becauase I'm not acutally sure what the units are, possibly pixels
    # lets work where if there was only 1 sequence it would be a size of 1000
    # therefore each vertice will be a node size that is the re_abund * 1000
    # NB the size does appear to be something similar to pixels
    vertices_sizes = [int(vertice_id_to_rel_abund_dict[vert] * 1000) for vert in vertices]

    # the edge list should be a list of tuples
    # currently we are not taking into account the length of the vertices
    edges_list = []
    for i in range(len(network_block)):
        if network_block[i] == 'EDGES':
            # for each line in the translate section
            for j in range(i + 1, len(network_block)):
                if network_block[j] == ';':
                    # then we have reached the end of the translation section
                    break
                items = network_block[j].replace(',', '').split(' ')[1:]
                edges_list.append((items[0], items[1]))

    # THis is useful, scroll down to where the 'def draw_networkx(......' is
    # https://networkx.github.io/documentation/networkx-1.9/_modules/networkx/drawing/nx_pylab.html#draw_networkx

    # so some of the vertice size are coming back as 0 because they were found in such low abundances.
    # because we are working with a subsampling of 1000 we should probably get rid of any of the vertices that
    # have size 0. However, we can't do that at this point because then we will be breaking the network
    # as such we are going to have to do this before we calculate the networks. This could actually work in our favour
    # as it could make the network less computationally expensive.

    g.add_nodes_from(vertices)
    g.add_edges_from(edges_list)
    # we should be able to
    # f, ax = plt.subplots(1, 1)
    # we will need to set the limits dynamically as we don't know what the positions are going to be.
    # I think we will want to have the x and y limits the same so that we end up with a square
    # we will therefore just look for the bigest and smallest values, add a buffer and then set the
    # axes limits to thsee
    spring_pos = nx.spring_layout(g)
    max_ax_val = 0
    min_ax_val = 9999
    for vert_pos_array in spring_pos.values():
        for ind in vert_pos_array:
            if ind > max_ax_val:
                max_ax_val = ind
            if ind < min_ax_val:
                min_ax_val = ind

    ax.set_ylim(min_ax_val - 0.2, max_ax_val + 0.2)
    ax.set_xlim(min_ax_val - 0.2, max_ax_val + 0.2)
    # to set the edge colour of the nodes we need to draw them seperately
    # nx.draw_networkx_nodes(g, pos=spring_pos, node_size=vertices_sizes, with_labels=False, alpha=0.5, edgecolors='black')
    nx.draw_networkx(g, pos=spring_pos, arrows=False, ax=ax, node_color=colour_for_network, alpha=1.0, node_size=vertices_sizes, with_labels=False )

    # https://stackoverflow.com/questions/22716161/how-can-one-modify-the-outline-color-of-a-node-in-networkx
    #https://matplotlib.org/api/collections_api.html
    ax.collections[0].set_edgecolor('grey')
    ax.collections[0].set_linewidth(0.5)
    ax.collections[1].set_linewidth(0.5)
    ax.set_axis_off()
    ax.set_title(type_id)
    apples = 'asdf'

def draw_network_split_colours(splits_tree_out_file, count_id_to_abund_dict, ax, colour_for_network, type_id):
    network_colour_dict = pickle.load(open('network_colour_dict.pickle', 'rb'))

    # networkx graph object
    g = nx.Graph()



    # we can work with a minimised version of the splits out file that is only the network block
    for i in range(len(splits_tree_out_file)):
        if splits_tree_out_file[i] == 'BEGIN Network;':
            network_block = splits_tree_out_file[i:]
            break

    # here we can now work with the network_block rather than the larger splits_tree_out_file
    # get a list of the nodes
    # these can be got from the splits_tree_outfile
    # we should use the numbered system that the splits tree is using. We can get a list of the different sequences
    # that each of the taxa represent form the translate part of the network block.
    # This way we can use this to run the sequences against SP to get the names of the sequencs
    # we can also make the 'translation dictionary' at the same time
    vertice_id_to_seq_id_dict = defaultdict(list)
    for i in range(len(network_block)):
        if network_block[i] == 'TRANSLATE':
            # for each line in the translate section
            for j in range(i + 1, len(network_block)):
                if network_block[j] == ';':
                    # then we have reached the end of the translation section
                    break
                items = network_block[j].replace('\'', '').replace(',', '').split(' ')
                vertice_id_to_seq_id_dict[items[0]].extend(items[1:])

    vertices = list(vertice_id_to_seq_id_dict.keys())

    # here we will use the dictionary that was created in the previous method (passed into this one)
    # to link the count id to the actual sequence, which can then be used to look up the abundance
    # This way we can have a vertice_id_to_rel_abund dict that we can use to put a size to the nodes

    vertice_id_to_rel_abund_dict = {}
    for vert_id in vertices:
        count_id_list = vertice_id_to_seq_id_dict[vert_id]
        if len(count_id_list) == 1:
            rel_abund = count_id_to_abund_dict[count_id_list[0]]
            vertice_id_to_rel_abund_dict[vert_id] = rel_abund
        elif len(count_id_list) > 1:
            # then we need sum the abundances for each of them
            tot = 0
            for count_id in count_id_list:
                rel_abund = count_id_to_abund_dict[count_id]
                tot += rel_abund
            vertice_id_to_rel_abund_dict[vert_id] = tot


    # the sizes are a little tricky to work out becauase I'm not acutally sure what the units are, possibly pixels
    # lets work where if there was only 1 sequence it would be a size of 1000
    # therefore each vertice will be a node size that is the re_abund * 1000
    # NB the size does appear to be something similar to pixels
    vertices_sizes = [int(vertice_id_to_rel_abund_dict[vert] * 1000) for vert in vertices]

    # the edge list should be a list of tuples
    # currently we are not taking into account the length of the vertices
    edges_list = []
    for i in range(len(network_block)):
        if network_block[i] == 'EDGES':
            # for each line in the translate section
            for j in range(i + 1, len(network_block)):
                if network_block[j] == ';':
                    # then we have reached the end of the translation section
                    break
                items = network_block[j].replace(',', '').split(' ')[1:]
                edges_list.append((items[0], items[1]))

    # THis is useful, scroll down to where the 'def draw_networkx(......' is
    # https://networkx.github.io/documentation/networkx-1.9/_modules/networkx/drawing/nx_pylab.html#draw_networkx

    # so some of the vertice size are coming back as 0 because they were found in such low abundances.
    # because we are working with a subsampling of 1000 we should probably get rid of any of the vertices that
    # have size 0. However, we can't do that at this point because then we will be breaking the network
    # as such we are going to have to do this before we calculate the networks. This could actually work in our favour
    # as it could make the network less computationally expensive.

    g.add_nodes_from(vertices)
    g.add_edges_from(edges_list)


    # we will need to set the limits dynamically as we don't know what the positions are going to be.
    # I think we will want to have the x and y limits the same so that we end up with a square
    # we will therefore just look for the bigest and smallest values, add a buffer and then set the
    # axes limits to thsee
    spring_pos = nx.spring_layout(g)
    max_ax_val = 0
    min_ax_val = 9999
    for vert_pos_array in spring_pos.values():
        for ind in vert_pos_array:
            if ind > max_ax_val:
                max_ax_val = ind
            if ind < min_ax_val:
                min_ax_val = ind

    ax.set_ylim(min_ax_val - 0.2, max_ax_val + 0.2)
    ax.set_xlim(min_ax_val - 0.2, max_ax_val + 0.2)
    # to set the edge colour of the nodes we need to draw them seperately
    # nx.draw_networkx_nodes(g, pos=spring_pos, node_size=vertices_sizes, with_labels=False, alpha=0.5, edgecolors='black')
    nx.draw_networkx(g, pos=spring_pos, arrows=False, ax=ax, node_color=colour_for_network, alpha=1.0, node_size=vertices_sizes, with_labels=False )

    # https://stackoverflow.com/questions/22716161/how-can-one-modify-the-outline-color-of-a-node-in-networkx
    #https://matplotlib.org/api/collections_api.html
    ax.collections[0].set_edgecolor('grey')
    ax.collections[0].set_linewidth(0.5)
    ax.collections[1].set_linewidth(0.5)
    ax.set_axis_off()
    ax.set_title(type_id)
    apples = 'asdf'

model_its_type_profiles()

