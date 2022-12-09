
# ------------------------------------
#
# To execute a command, place the cursor on that command's line
# and hit Ctrl + Enter. Try executing the following command:

print("I am executing this")

# You should see that, in the console, the message "I am executing this"
# was printed. In general, the results of any command you execute
# will appear in the console.
#
# ------------------------------------
# 
# Let us import the artifactor package into this project by executing

library(artifactor)

# You  will be asked whether to install Python dependencies.
# Write "Yes" in the terminal to accept this.  
# 
# ---------------DATA LOADING---------------------
# 
#
#
# Now we can start using our scientific EEG package. First of all, let's
# read the EEG data you've exported from EDFReader. 
# To load EEG data, execute the following command.

eeg <- load_eeg("test_eeg.csv")

# Now the variable "eeg" holds an eeg object. You can access its data by
# writing "eeg@data". For example, execute the following command
# and see what you get.
# 

eeg@data

# You should see a database appear in the console. The first column is Time, 
# and the rest are EEG channels. This database contains the EEG measurements
# over a night of sleep from a subject.
# 
# ---------------DATA CLEANING---------------------
# 
# You should also see some values in the database are red and read "NA". 
# These are missing values, and analysis can not be conducted with them.
# For the purpose of this test, we also want to remove the Chin, ECG and Leg
# channels, so we are left only with the information gathered from the 
# electrodes in the scalp. So go ahead and execute the following commands.
# 

eeg@data[10:14] <- NULL # This erases columns 10, 11, 12, 13, 14.
                        # These are the columns corresponding to non-cranial electrodes.

eeg@data <- na.omit(eeg@data) # This removes NA values.
                        
# Now inspect the data once more. You should see no more NA values.
# 

eeg@data
 
# ---------------RESAMPLING AND PARTITIONING THE EEG---------------------
# 
# If ou look at the eeg@data, you'll see we have one measurement every 0.002
# seconds. This might be an unnecesarily high resolution for our analysis,
# so we will resample the EEG. Resampling means keeping only one every n values,
# where n is a natural number. Execute the following command
# 

resampled_eeg <- resample_eeg(eeg, 5)

# Now we have a new EEG object called resampled_eeg, and if you plot its data
# you'll see it has one measurement every 0.01 seconds.

resampled_eeg@data


# ---------------SEEING THE EEG---------------------
# 
# We can draw our EEG at specific time intervals to get a sense of what
# we area dealing with. There are two types of plots: interactive and static
# plots. We will call both of them.
#
# Say we wanted to see the fifth minute of our EEG record. To do this,
# we need to take the subset of our EEG in that time, and then plotting 
# that subset. We can do this in one line.
 
fifth_minute <- subset_eeg(resampled_eeg, 300, 360)
plot(fifth_minute)

# You should see a plot appear in the RStudio "Plots" section.
# Click "Zoom", right above the plot, to take a better look.
# 
# The function `subset_eeg` takes an eeg and two numbers representing
# time in seconds. It returns the section of that EEG that's between the
# first and second time in seconds. The function `plot_eeg` simply plots an 
# eeg object. 
# 
# So, what we did before was
# a) store the subset of the resampled eeg from time 300 to 360 into a variable 
# called `fifth_minute`.
# b) plot the eeg that we stored into that variable.
# 
# Using the same idea, go ahead and try to plot the EEG from the tenth
# to the eleventh minute. That would from seconds 600 to 660.
# To do this, fill in the blanks! And use the code above to help yourself.

tenth_minute <- ?
plot_eeg(?)


# ---------------ANALYZING THE EEG ---------------------
# 
# -------------------------------------------------------

# Execute this:

rs <- create_epoch_data()
sbs <- subset_eeg(resampled_eeg, 0, 600)
an <- artf_stepwise(sbs, 60)

# Make a static plot of the analysis:
# does it appear in the Plots pannel?
plot(an)

# The following command may take a while to complete.
# After completion, a web page showing an interactive plot
# of the analysis should open in your default browser.
# If you can, take a screenshot of what it looks like, save it
# and send it to me.
iplot_analysis(set_for_iplot(an), show=TRUE, save=FALSE)

# Execute this:
rs <- rs %>% update_epochs(an)
comparative_rs <- create_epoch_data()
comparative_an <- stepwise_analysis(sbs, 60, thresh=0, alpha=20)
comparative_rs <- comparative_rs %>% update_epochs(comparative_an)

# The following commands will save two files into your project's directory,
# named "rs.csv" and "crs.csv". Send those files to me once they're generated
# so I can compare your results (they should be the same).

write_csv(rs, "rs.csv")
write_csv(comparative_rs, "crs.csv")
