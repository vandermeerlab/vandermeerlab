Shortcut session-specific information
=====================================

Each recording session contains a unique *_info.py file that has the following information:

* `session_id` : str
* `species` : str
* `behavior` : str
* `target` : list of str
* `experimenter` : str
* `task_times` : dict  
	Where prerecord, phase1, pauseA, phase2, pauseB, phase3, postrecord are keys and the value is a list of floats (start, stop times.)
* `pxl_to_cm` : tuple of floats (x, y conversion factor)
* `fs` : int or float
* `good_lfp` : list of str
* `good_swr` : list of str
* `good_theta` : list of str
* `path_pts` : dict  
	Where the key is a str representing a point in 2D space, and the value is a list of 2 ints, representing x, y location in pixels (before the conversion to cm is applied).
* `u_trajectory` : list of lists containing x, y floats
* `shortcut_trajectory` : list of lists containing x, y floats
* `novel_trajectory` : list of lists containing x, y floats

