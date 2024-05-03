# RSD2_Analysis
Since there is no documentation, I will try to reconstruct what the Turin group does here.

## Infrastructure
Project seems to have a main file called `Analysis_Digitizer_wfm.cpp` that serves as an entrypoint for the analysis.

For each data file, there must be an entry in the *Input Card* `Input_Folder_Digitizer_wfm.txt`.

The `bash_files` folder seems to contain shellscripts that start the analysis and broadcast its current status and results over a discord bot.

## Input card format
There appear to be several different keywords in the `.txt` file that serve as input settings to the analysis:

 - `Run` seems to set the run number,
 - `StartingTrig` seems to set the relative starting point for the trigger (where is the pulse within the readout window?),
 - `MaxEventPerPoint` sets the number of events per laser injection point,
 - `SkipPos` skips specific laser injection points,
 - `Nchannel` sets the number of channels to analyze,
 - `MaxEvent` seems to set the maximum number of events the analysis will look at (no limit is -1),
 - `ana` triggers the actual analysis and gets the input file path as a parameter,
 - `nmedia` ???
 - `Max_number_channels` set the maximum number of available channels (?)

For some reason, most of the entries seem to have a `Run -1` at the end... I don't know why, but we will use it, too, until we know what it does.

## KIT Modifications
In the following, I will list the modifications that needed to be done.

### `Analysis_Digitizer_wfm.cpp`
