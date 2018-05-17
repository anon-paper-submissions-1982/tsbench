# first line: 814
@_memory.cache
def _compute_ucr_snrs():
    from _python.datasets import ucr
    dsets = ucr.allUCRDatasets()

    snrs = {}

    for dset in dsets:
        X = dset.X
        minval, maxval = np.min(X), np.max(X)
        spread = maxval - minval
        for nbits in (8, 16):
            dtype = {8: np.uint8, 16: np.uint16}[nbits]
            scale = (1 << nbits) / spread
            X_quant = np.array((X - minval) * scale, dtype=dtype)
            # X_quant = np.array((X - minval) * scale + .5, dtype=dtype)

            X_hat = (X_quant.astype(np.float64) / scale) + minval

            diffs = X - X_hat

            snr = np.var(X.ravel()) / np.var(diffs.ravel())
            snrs[(dset.name, nbits)] = snr

    return snrs
