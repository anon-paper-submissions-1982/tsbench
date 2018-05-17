# first line: 221
@_memory.cache
def _read_and_clean_ucr_results(no_preprocs, results_path=cfg.UCR_RESULTS_PATH):
    df = pd.read_csv(results_path)

    # df = df[df['Order'] == 'c']
    df['Filename'] = df['Filename'].apply(lambda s: os.path.basename(s).split('.')[0])
    df = df.sort_values(['Filename', 'Algorithm'])
    # if others_deltas:
    #     df = df[df['Deltas'] | df['Algorithm'].str.startswith('Sprintz')]
    # else:
    #     df = df[~df['Deltas']]
    if no_preprocs:
        df = df[df['Preprocs'] == cfg.Preproc.NONE]
        df = df[['Nbits', 'Algorithm', 'Filename', 'Order', 'Ratio']]
    else:
        df = df[['Nbits', 'Algorithm', 'Filename', 'Order', 'Preprocs', 'Ratio']]
    df = df[df['Algorithm'] != 'Memcpy']
    df['Ratio'] = 100. / df['Ratio']  # bench reports % of original size
    return df
