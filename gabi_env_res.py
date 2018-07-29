# Now we want to answer the second question which is more a graphical exploration
# first get a dictionary that links coral sample to mucus sample
coral_info_name_to_muc_info_name_dict = get_coral_to_mucus_dict()

# now we need to convert the coral and mucus names as they are found in the info table to the names
# as they are found int he database.
coral_name_to_muc_name_dict = get_db_named_coral_to_mucus_dict(coral_info_name_to_muc_info_name_dict, data_set_samples, sample)
# here we have a dict that links the coral sample name to the associated mucus sample name.

# now we can start to get the numbers in order and then plot
# locations where the count table outputs are found.
std_out_dir = '/home/humebc/phylogeneticSoftware/SymPortal_interim_250318/outputs/25'
env_out_dir = '/home/humebc/phylogeneticSoftware/SymPortal_interim_250318/outputs/environmental_type_screening/86'

# The data that will hold the coral data
coral_div_df = pd.read_csv('{}/25_87.DIVs.absolute.txt'.format(std_out_dir), sep='\t', lineterminator='\n')
column_names_to_drop = ['noName Clade A',       'noName Clade B',       'noName Clade C',       'noName Clade D',
                        'noName Clade E',       'noName Clade F',       'noName Clade G',       'noName Clade H',
                        'noName Clade I',       'raw_contigs',  'post_qc_absolute_seqs',        'post_qc_unique_seqs',
                        'post_taxa_id_absolute_symbiodinium_seqs',      'post_taxa_id_unique_symbiodinium_seqs',
                        'post_taxa_id_absolute_non_symbiodinium_seqs',  'post_taxa_id_unique_non_symbiodinium_seqs']
coral_div_df_no_meta_no_noName = coral_div_df.drop(column_names_to_drop, axis=1)

# Then the environmental data too.
total_div_df = pd.read_csv('{}/25_86.DIVs.env.absolute.total.txt'.format(env_out_dir), sep='\t', lineterminator='\n')
total_div_df_no_meta_no_noName = total_div_df.drop(column_names_to_drop, axis=1)




# We are going to want all of the sequences to be comparable between the plots.
#It will be too much to work with all of the sequenecs coloured.
# I'm thinking we should colour the DIVs but then the non-DIVs should be in black and white / shades of grey
# list_of_env_headers = list(total_div_df)[1:202]
# list_of_coral_headers = list(coral_div_df)[1:161]
list_of_env_headers = list(total_div_df_no_meta_no_noName)[1:]
list_of_coral_headers = list(coral_div_df_no_meta_no_noName)[1:]
# We can do a better job of getting an ordered list of the headers (i.e. DIVs)
# but for the time being let's just work through the list_of_env_headers followed by the coral headers
# first get a set of the environmental header
all_divs = list_of_env_headers + list_of_coral_headers
unique_list_of_all_divs = ordered_unique_list(all_divs)

#TODO here we want to identify which of the DIVs are found at an abundance > 1% in the samples
# we will have to go through each of the columns of the tables and see what the proportion of the DIV is
# for this we can use the proportion dataSets
env_rel_abund_df_to_rem = pd.read_csv('{}/25_86.DIVs.env.proportion.total.txt'.format(env_out_dir), sep='\t', lineterminator='\n')
env_rel_abund_df = env_rel_abund_df_to_rem.drop(column_names_to_drop, axis=1)
env_rel_abund_df.drop(env_rel_abund_df.index[len(env_rel_abund_df) - 1], inplace=True)
coral_rel_abund_df_to_rem = pd.read_csv('{}/25_87.DIVs.relative.txt'.format(std_out_dir), sep='\t', lineterminator='\n')
coral_rel_abund_df = coral_rel_abund_df_to_rem.drop(column_names_to_drop, axis=1)
coral_rel_abund_df.drop(coral_rel_abund_df.index[len(coral_rel_abund_df) - 1], inplace=True)

### N.B. this cutoff can be changed to speed up the plotting. i.e. for the time beint we can work with 5%
cutoff = 0.05

# get a list of the clades found.
# we can go through the env divs first
clades_found = set()
div_names_above_one_percent = []
for header in list(env_rel_abund_df)[1:]:
    if '_' in header:
        clades_found.add(header.split('_')[1])
    else:
        clades_found.add(header[0])
    if max(list(env_rel_abund_df[header])) >= cutoff:
        # Then this DIV is found at > 1 % in one of the samples
        div_names_above_one_percent.append(header)

# we can go through the coral samples now
for header in list(coral_rel_abund_df)[1:]:
    if '_' in header:
        clades_found.add(header.split('_')[1])
    else:
        clades_found.add(header[0])
    if header not in div_names_above_one_percent:
        if max(list(coral_rel_abund_df[header])) >= cutoff:
            div_names_above_one_percent.append(header)

clades_found = sorted(clades_found)
# Here we have a list of DIVs that were found above 1% in at least one of the samples
# hopefully we have lowered our DIV number to lower than 4000 divs now.

#TODO now we need to do the clade grouping of the sequences that were below the 1% mark.
# Do this on a sample by sample basis. Where we start to add up the div abundances.
# for each of the samples have an additional set of clades


# This is the list of 269 colours in hex
colour_list_template = get_colour_list()
grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
# From this we need to get the actual order of colours we want.
# we will colour the DIVs in colour and the noName sequences in black and white greys.
colours_in_order_of_DIVs = []
colour_index = 3
b_w_index = 0
for DIV in div_names_above_one_percent:
    # if '_' in DIV:
    #     # then this is a noName seq and it should be one of the Greys
    #     colours_in_order_of_DIVs.append(grey_palette[b_w_index%len(grey_palette)])
    #     b_w_index += 1
    # else:
    #     # then this is a named seq and it should be in colour
    #     colours_in_order_of_DIVs.append(colour_list_template[colour_index%len(colour_list_template)])
    #     colour_index += 1
    # here we are just going to work on colouring all seqs but the clade groupings will be grey colours
    colours_in_order_of_DIVs.append(colour_list_template[colour_index%len(colour_list_template)])
    colour_index += 1


# now add the clade groupings colours
for clade_found in clades_found:
    colours_in_order_of_DIVs.append(colour_list_template[colour_index%len(colour_list_template)])
    colour_index += 1
f, axarr = plt.subplots(4, 1, sharey='row', sharex='col')
axarr_index = -1

# go sample type by sample type and have a subplot for each
for smpl_type in ['turf', 'coral', 'mucus', 'seawater']:
    # find which count table we should be using for the look up
    if smpl_type == 'coral':
        df_to_use = coral_div_df
    else:
        df_to_use = total_div_df

    # the holder for the DIV abundances for each sample which will have its own list
    data_for_plot_for_transposing = []
    samples_to_plt = []

    # get a list of the samples that we're going to be working with for this sample type
    for k, v in sample_name_to_sample_type_dict.items():
        if v == smpl_type:
            samples_to_plt.append(k)
    axarr_index += 1
    # here we have a list of samples of the chosen sample type

    # TODO now we need to order the samples according to similarity


    # use the count table to look up the abundances of each of the divs above the 1% cutoff
    for smpl_name in samples_to_plt:
        print('{}:{}'.format(smpl_type, smpl_name))
        # get which row we're working with i.e. which sample so that we can grab the seq abund data
        smpl_index = list(df_to_use['Samples']).index(smpl_name)
        # add data for each of the DIVs

        # at first we will grab the abosolute abundances and then convert to rel abundance after
        sample_data_to_norm = []
        clade_grouping_count_to_norm = {clade: 0 for clade in clades_found}
        # for every DIV check to see if it is one of the  1%.
        # If so add the abundance if possible, or 0
        # else, add this abundance to the clade groupings.
        for DIV in unique_list_of_all_divs:
            if DIV in div_names_above_one_percent:
                try:
                    sample_data_to_norm.append(df_to_use[DIV][smpl_index])
                except:
                    sample_data_to_norm.append(0)
            else:
                if '_' in DIV:
                    # then this is a noName
                    of_clade = DIV.split('_')[1]
                else:
                    # then this is a names
                    of_clade = DIV[0]
                try:
                    clade_grouping_count_to_norm[of_clade] += df_to_use[DIV][smpl_index]
                except:
                    clade_grouping_count_to_norm[of_clade] += 0
        tot = sum(sample_data_to_norm) + sum(clade_grouping_count_to_norm.values())

        smpl_total_data = [a/tot for a in sample_data_to_norm]
        clade_grouping_count = [clade_grouping_count_to_norm[clade]/tot for clade in clades_found]
        # now add both the names div data and the clade groupings to the dataholder.
        # append the clade groupings to the div counts in the order of the clades list
        data_for_plot_for_transposing.append(smpl_total_data + clade_grouping_count)

    # now we have the relative DIV abundances for each of the samples. A list of each sample.
    # we must now transpose this list of lists so that we have a list for each of the divs to be plotted.
    data_for_plot = list(map(list, zip(*data_for_plot_for_transposing)))

    # # for seq in sortedSeqs:
    # #         barList.append(ax.bar(ind, plotInfoDict[seq], width, bottom=bottom, color=c_palette[colourCounter % len(c_palette)]))
    # #         bottom = [L + M for L, M in zip(bottom, plotInfoDict[seq])]
    # #         colourCounter += 1
    # #     ax.legend(handles=barList, labels=[seq for seq in sortedSeqs], loc='center left', bbox_to_anchor=(1, 0.5))
    # #     ax.set_title(str(at))
    # # plt.tight_layout()
    # # plt.savefig('debug.svg', format="svg")



    # TODO this is way too many sequences, we need to cut this down, I think we're going to have to group those
    # sequences grouped below 1% on average across all samples (or maybe within sample) into clade.
    # We also need to order the samples according to similarity.
    # https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html
    # The above looks great for this.

    bottom = [0 for k in range(len(samples_to_plt))]
    bar_list = []
    len_div = len(div_names_above_one_percent) + len(clades_found)
    for i in range(len_div):
        if i < len(div_names_above_one_percent):
            print('DIV:{}:{}/{}'.format(unique_list_of_all_divs[i], i, len_div))
        else:
            # print('DIV:clade {} grouping :{}/{}'.format(clades_found[i-len(unique_list_of_all_divs)], i, len_div))
            print('Printing clade groups')
        bar_list.append(axarr[axarr_index].bar(range(len(samples_to_plt)), data_for_plot[i], 1, bottom, color=colours_in_order_of_DIVs[i]))
        bottom = [L + M for L, M in zip(bottom, data_for_plot[i])]
f.savefig('/home/humebc/projects/env_res/fig1_seq_diversity.svg')
f.savefig('/home/humebc/projects/env_res/fig1_seq_diversity.png')


# totals
apples = 'asdf'

apples = 'asdf'
return
