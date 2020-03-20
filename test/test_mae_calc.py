import numpy as np
from wrf_ens_tools import calc

# Initialize test data
predictions = np.ones((50, 50)) * 3. # forecast
targets = np.ones((50, 50)) * 1.     # observations
expected_mae = 2.

def test_mae_default():
    # Test mean absolute error calculation assuming optional arg defaults
    mae = calc.mae(predictions=predictions, targets=targets,
                    axis=None, nan=False)
    assert(mae == expected_mae)
