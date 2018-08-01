import pandas as pd
import sys
import pickle
import os
import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import statistics


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
            top_n_seqs_done = []
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

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='s', s=8, c='r')
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
                    axarr[axarr_index].set_xlabel('pre-MED / post-MED')
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

                ax2.scatter(x=ind[ind_index] + 0.375, y=mean, marker='s', s=8, c='r')
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

def make_venns():
    ''' Here we will aim to make venn diagrams we will make a set of three 2 cirle venns which will
    just show the comparison of the env_samples to the coral.
    We will also make a set of four 3 circle venn diagrams which will show all 3 way combinations.
    There is a neat little module called matplotlib_venn which we can use to do this and it is dead simple.
    All you need to do is pass it a list of sets.'''

    # read in the dataframes created previously
    sp_output_df = pickle.load(open('sp_output_df.pickle', 'rb'))
    QC_info_df = pickle.load(open('QC_info_df.pickle', 'rb'))
    info_df = pickle.load(open('info_df.pickle', 'rb'))

    # create a dictionary that is sample name to env_type
    sample_name_to_sample_type_dict = {}
    for info_index in info_df.index.values.tolist():
        if info_df.loc[info_index, 'environ type'] == 'coral':
            sample_name_to_sample_type_dict[info_df.loc[info_index, 'sample_name']] = 'coral'
        elif info_df.loc[info_index, 'environ type'] == 'sea_water':
            sample_name_to_sample_type_dict[info_df.loc[info_index, 'sample_name']] = 'sea_water'
        elif info_df.loc[info_index, 'environ type'] == 'sed_far':
            sample_name_to_sample_type_dict[info_df.loc[info_index, 'sample_name']] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'sed_close':
            sample_name_to_sample_type_dict[info_df.loc[info_index, 'sample_name']] = 'sed'
        elif info_df.loc[info_index, 'environ type'] == 'mucus':
            sample_name_to_sample_type_dict[info_df.loc[info_index, 'sample_name']] = 'mucus'
        elif info_df.loc[info_index, 'environ type'] == 'turf':
            sample_name_to_sample_type_dict[info_df.loc[info_index, 'sample_name']] = 'turf'

    # here we have the dict populated and we can now get to work pulling out the set of sequence
    # names found in each of the sample types
    sample_types = ['coral', 'sea_water', 'sed', 'mucus', 'turf']

    set_list = [set() for type in sample_types]

    # now work through the sp output and check to see which columns are non-zero columns and add these column
    # labels into the respective sets
    for sp_output_index in sp_output_df.index.values.tolist():
        type_of_sample  = sample_name_to_sample_type_dict[sp_output_df.loc[sp_output_index, 'sample_name']]
        non_zero_seqs = list(sp_output_df.loc[sp_output_index][sp_output_df.loc[sp_output_index] > 0].index)
        set_list[sample_types.index(type_of_sample)].update(non_zero_seqs)

    # here we should have the sets populated
    # now create a dictionary that is the sample_type name to the set

    # then simply do permutations for the 4 x 3

    # and maybe manually do the three venns using the sets
    



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

def get_db_named_coral_to_mucus_dict(coral_info_name_to_muc_info_name_dict, data_set_samples, sample):
    coral_name_to_muc_name_dict = {}
    for k, v in coral_info_name_to_muc_info_name_dict.items():
        coral_count = 0
        coral_name_match = str()
        mucus_count = 0
        mucus_name_match = str()
        for sample in data_set_samples:
            # look up in the info table what the sample type is and append it to a dict
            if sample.name.replace('-', '_') in k or k in sample.name.replace('-', '_'):
                # Then we have a match for the info name and the actual sample name for the coral
                coral_count += 1
                coral_name_match = sample.name
            if sample.name.replace('-', '_') in v or v in sample.name.replace('-', '_'):
                # Then we have a match for the info name and the actual sample name for the mucus
                mucus_count += 1
                mucus_name_match = sample.name

        if coral_count == 1 and mucus_count == 1:
            # then we only found one match and we are good to continue.
            # get the sample_type for this sample and add this to the sample_name_to_sample_type_dict
            coral_name_to_muc_name_dict[coral_name_match] = mucus_name_match
        else:
            sys.exit('Count != 1 {}'.format(sample.name))
    return coral_name_to_muc_name_dict

def get_coral_to_mucus_dict():
    coral_to_mucus_info_df = pd.read_csv('/home/humebc/projects/env_res/seq_info_coral_mucus_link.csv', sep=',',
                                         lineterminator='\n')
    coral_info_name_to_muc_info_name_dict = {}
    for i in range(1, 69):
        mucus_found = False
        coral_found = False
        coral_sample_name = ''
        mucus_sample_name = ''
        for j, df_name in enumerate(coral_to_mucus_info_df['Sample no.']):

            if df_name.strip() == 'coral {}'.format(i):
                # then this is the coral sample we are intersted in
                coral_sample_name = coral_to_mucus_info_df['Sample Name'][j]
                coral_found = True
            if df_name.strip() == 'mucus {}'.format(i):
                # then this is the paired mucus sample
                mucus_sample_name = coral_to_mucus_info_df['Sample Name'][j]
                mucus_found = True

        if coral_found and mucus_found:
            # Then we can add this to the dict
            coral_info_name_to_muc_info_name_dict[coral_sample_name] = mucus_sample_name
        else:
            # then we just wont work with this pair.
            print('{} not found'.format(i))
            continue

    return coral_info_name_to_muc_info_name_dict

make_venns()