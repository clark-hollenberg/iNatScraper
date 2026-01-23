# iNaturalist API tool for the NatureServe Network
Written by Clark Hollenberg in 2025, adapting code from Kyle Kaskie at MTNHP.  
**Basic Python knowledge is currently necessary to configure and run this tool. No major code edits should be required.**

## What it does:
* Queries iNaturalist for observations matching a user defined tracking list within your jurisdiction.
* Provides access to obscured coordinates for users which have trusted your iNaturalist project.
* Filters observations based on location accuracy and whether or not they have an expert ID.

## How to configure it:
### If querying for private coordinates observations based on project trust
* iNaturalist project with trusting users
* iNaturalist app creation (by owner of the project):
    * Meet Activity Thresholds: You must have 10+ "improving" identifications (IDs that change a taxon from a higher level to a lower level, e.g., Insecta to Species) on other users' observations in the past 30 days.
    * Account Age: Your account must be older than 60 days.
    * Access the Application Form: Go to iNaturalist.org/oauth/applications and click "New Application".
    * input your password information input into the config.json (make a copy)
### Required inputs:
* Tracking list which contains scientific names and ELCODEs. Note that you should filter the tracking list before each run to what group you are interested in. It is possible to do the entire tracking list at once but make sure this is necessary.
* List of iNaturalist experts (provided in GitHub repo. Modified from list from MTNHP.)
* State code (number) in iNaturalist
    * found by searching your state on the explore page and viewing URL. Ex: Colorado = 34. https://www.inaturalist.org/observations?place_id=34.
### Optional input:
* Date last searched for each species. We use this as a way to track which records we've already reviewed for Biotics. There are likely cleaner ways to do this.

## How to run it:
* git clone 
* configure Python interpreter
* copy and update config.json file with your info
* run build_taxon_ids.py
  * review output for any var. or ssp. input names
  * recommended to review output for species complexes or other issues (see iNat names column)
* run CNHP_iNat_Scraper.py
* view saved files in tracked_obs folder
  * observations_for_review contains all obs with unobscured coordinates
    * note that there are expert review columns that can be filtered on to see which have been ID'd by an expert.
  * observations_needing_accuracy contains all other tracked obs which are obscured (consider reaching out to iNat observer for trust in project)
  


