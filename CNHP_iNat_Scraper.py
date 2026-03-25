"""
iNaturalist API download for tracked species

Created by Clark Hollenberg in 2025.

Downloads iNaturalist observations in Colorado that are found on the supplied tracking list.
Identifies which observations contain high-res location data.
Identifies which observations have been expert identified.

Requires configuration in config.json file and the creation of an iNaturalist API "app".
Users must "trust" the project in iNaturalist for you to view private coordinates.
The tracking list can be downloaded from Biotics to refresh the desired query, but 
the taxon_id in iNaturalist must be identified for each scientific name, which requires API
queries. The iNaturalist system usually handles synonymy correctly, but there are some species
which may need manual resolution. For this reason, the taxon_name_mappings must be handled carefully
to preserve manual user inputs and not overwrite taxon_ids from iNaturalist which are incorrect.

Note that we use a "date_checked.csv" file based on "Collections searched" in Biotics, where we record
when all the records in iNat were "caught up" this is used as a created_after date filter in the query if applicable.

"""
# .venv\Scripts\activate
import requests
import pandas as pd
import numpy as np
import json
import os
import datetime
import time
from typing import List, Dict, Tuple

class iNatScraper:
    def __init__(self, config_file="config_private.json"):
        """Initialize with configuration file"""
        self.config_file = config_file
        self.config = self.load_config()
        self.tracking_list = self.load_tracking_list()
        self.taxon_name_map = self.load_taxon_name_map() 
        self.taxon_date_map = self.load_date_checked()
        self.timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        
    def load_config(self) -> Dict:
        with open(self.config_file, 'r') as f:
            return json.load(f)
    
    def load_experts(self) -> pd.DataFrame:
        return pd.read_csv(self.config['experts']['csv_file'], encoding='latin1')
    
    def load_tracking_list(self) -> pd.DataFrame:

        ### here is where tracking list can be filtered.
        ### any new names must be compiled in build_taxon_ids.py
        tracking_list = pd.read_csv(self.config['taxonomy']['tracking_list'], encoding='latin1')
        return tracking_list
    
    def load_taxon_name_map(self) -> pd.DataFrame:
        taxon_name_map = pd.read_csv(self.config['taxonomy']['taxon_name_map'])
        taxon_name_map = taxon_name_map[taxon_name_map["taxon_id"].notna()]
        return taxon_name_map

    def load_date_checked(self) -> Dict[str, str]:
        """
        Load date_checked.csv and create mapping from taxon_id to last checked date
    
        """
        date_checked_file = self.config["taxonomy"]["date_checked"]
        
        if not os.path.exists(date_checked_file) or date_checked_file=="":
            print(f"Date checked file not found: {date_checked_file}")
            print("All observations will be downloaded from beginning of time.")
            return {}
        
        try:
            date_df = pd.read_csv(date_checked_file)
            print(f"Loaded date_checked.csv with {len(date_df)} entries")
            return date_df
            
        except Exception as e:
            print(f"Error loading date_checked.csv: {e}")
            return {}

    def get_oauth_access_token(self) -> str:
        """Get OAuth2 access token for private coordinate access"""
        auth_config = self.config.get('authentication', {})
        username = auth_config.get('username', '').strip()
        password = auth_config.get('password', '').strip()
        app_id = auth_config.get('app_id', '').strip()
        app_secret = auth_config.get('app_secret', '').strip()
        
        if not all([username, password, app_id, app_secret]):
            print("Error: Missing OAuth2 credentials in config file")
            print("Required: username, password, app_id, app_secret")
            return ""
        
        try:
            site = "https://www.inaturalist.org"
            payload = {
                'client_id': app_id,
                'client_secret': app_secret,
                'grant_type': "password",
                'username': username,
                'password': password
            }

            # Get OAuth access token
            access_token = requests.post(f"{site}/oauth/token", data=payload).json()["access_token"]
            headers = {"Authorization": f"Bearer {access_token}"}
            # Get API token
            api_token = requests.get(f"{site}/users/api_token", headers=headers).json()["api_token"]
            return api_token

        except Exception as e:
            print(f"Error getting OAuth2 access token: {e}")
            print("Check your credentials in config.json")
            return ""
    
    

    def chunks(self, lst: List, n: int):
        """Helper function to yield successive n-sized chunks from lst"""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def download_observations(self) -> pd.DataFrame:
        """Download iNaturalist observations with private coordinates using v1 API + Bearer token"""
        print("Downloading iNaturalist observations with private coordinates (v1 API + Bearer token)...")

        # Get OAuth2 access token (JWT string)
        api_token = self.get_oauth_access_token()
        if not api_token:
            print("Error: Could not obtain OAuth2 access token")
            return pd.DataFrame()

        headers = {"Authorization": f"Bearer {api_token}"}

        # Filter taxon mappings to the supplied tracking list
        # filter to species with taxon_id
        # Filter valid taxon IDs and only those in tracking list
        taxon_id_map = self.taxon_name_map[self.taxon_name_map["taxon_id"].notna()]
        taxon_id_map = taxon_id_map[taxon_id_map["S_ELMT_ID"].isin(self.tracking_list["S_ELMT_ID"])]
        taxon_id_map.to_csv(f"tracked_obs/SearchedSpecies_{self.timestamp}.csv", index=False) # save list of species searched

        # Collect list of scientific names being tracked
        scientific_names = taxon_id_map["SNAME"].unique().tolist()
        print(f"Tracking {len(scientific_names)} species")

        # Precompute S_ELMT_ID lookup set for faster matching
        if self.taxon_date_map != {}:
            date_id_set = set(self.taxon_date_map["S_ELMT_ID"].astype(str))
        else:
            date_id_set = []

        # Initialize output containers
        individual_requests = []
        no_date_filter_species = []
        missing_species = []

        # Iterate through rows to build requests
        for _, row in taxon_id_map.iterrows():
            taxon_id_str = str(int(row["taxon_id"]))
            scientific_name = row["SNAME"]
            element_ID = str(int(row["S_ELMT_ID"]))

            if element_ID in date_id_set:
                # Retrieve the matching date_checked value
                date_filtered = self.taxon_date_map.loc[
                    self.taxon_date_map["S_ELMT_ID"].astype(str) == element_ID, "date_checked"
                ]
                if date_filtered.empty:
                    missing_species.append((taxon_id_str, scientific_name))
                    continue

                inat_date = date_filtered.iloc[0]

                # Try to merge with an existing request (same date)
                merged = False
                for idx, (existing_ids, existing_date, existing_names) in enumerate(individual_requests):
                    if existing_date == inat_date:
                        # check to make sure that the id batch is not getting too long
                        if len(str(existing_ids).split(",")) >= 50:
                            #continue to create new id batch
                            continue

                        existing_id_set = {s.strip() for s in str(existing_ids).split(",") if s.strip()}
                        existing_name_set = {s.strip() for s in str(existing_names).split(";") if s.strip()}

                        # Append only if not already present
                        if taxon_id_str not in existing_id_set:
                            existing_ids = ",".join(sorted(existing_id_set | {taxon_id_str}))
                        if scientific_name not in existing_name_set:
                            existing_names = ";".join(sorted(existing_name_set | {scientific_name}))

                        individual_requests[idx] = (existing_ids, existing_date, existing_names)
                        merged = True
                        break
                            

                if not merged:
                    individual_requests.append((taxon_id_str, inat_date, scientific_name))
            else:
                no_date_filter_species.append((taxon_id_str, scientific_name))

        # Warn if any species weren’t found in date map
        if missing_species:
            print(f"Warning: {len(missing_species)} species not found in date_map")

        all_observations = []

        # Helper to fetch pages for given params
        def fetch_pages(params, query_flag):
            page = 1
            per_page = params.get('per_page', 200)
            base_params = dict(params)
            while True:
                params['page'] = page
                try:            
                    r = requests.get('https://api.inaturalist.org/v1/observations', headers=headers, params=params, timeout=30)
                    r.raise_for_status()
                    data = r.json()
                    results = data.get('results', [])

                    if not results:
                        break
                    
                    # flags for debugging
                    for obs in results:
                        obs["_query_taxon_id"] = base_params.get("taxon_id")
                        obs["_query_flag"] = query_flag
                        obs["_query_updated_since"] = base_params.get("updated_since")
                        obs["_query_place_id"] = base_params.get("place_id")
                        obs["_query_quality_grade"] = base_params.get("quality_grade")

                    all_observations.extend(results)
                    time.sleep(0.5)

                    if len(results) < per_page:
                        break
                    page += 1
                    time.sleep(0.5)

                except Exception as e:
                    print(f"Error fetching page. Error: {e}")

        # Individual requests (with created_after date filters)
        for taxon_id, start_date, scientific_name in individual_requests:
            print(f"Processing {scientific_name} (ID: {taxon_id}) with date filter...")
            # Parse start_date if present
            created_after = None
            if start_date:
                for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%m-%d-%Y"):
                    try:
                        created_after = datetime.datetime.strptime(start_date, fmt).date()
                        break
                    except ValueError:
                        continue
                if created_after is None:
                    print(f"  Warning: Could not parse date '{start_date}' for {scientific_name}")
            params = {
                'taxon_id': str(taxon_id),
                'place_id': self.config['search_area']['place_id'],
                'project_id': self.config["authentication"]["project"],
                'quality_grade': ','.join(self.config['search_area'].get('quality_grade', [])) if isinstance(self.config['search_area'].get('quality_grade'), list) else self.config['search_area'].get('quality_grade'),
                'per_page': self.config['search_area'].get('per_page', 200),
                'order_by': 'created_at',
                'order': 'desc'
            }
            if created_after:
                params['updated_since'] = created_after.isoformat()

            fetch_pages(params, "date_filtered")
            time.sleep(0.5)

        # Batched requests for species without date filters (use comma-separated taxon_id)
        batch_size = 50
        batches = list(self.chunks(no_date_filter_species, batch_size))
        for batch_num, batch in enumerate(batches, start=1):
            taxon_ids = [str(tid) for tid, _ in batch]
            print(f"Processing batch {batch_num}/{len(batches)}: {len(taxon_ids)} species...")
            params = {
                'taxon_id': ','.join(taxon_ids),
                'place_id': self.config['search_area']['place_id'],
                'project_id': self.config["authentication"]["project"],
                'quality_grade': ','.join(self.config['search_area'].get('quality_grade', [])) if isinstance(self.config['search_area'].get('quality_grade'), list) else self.config['search_area'].get('quality_grade'),
                'per_page': self.config['search_area'].get('per_page', 200),
                'order_by': 'created_at',
                'order': 'desc'
            }

            fetch_pages(params, "batch_filtered")
            time.sleep(1)  # be nice to API between batches

        print(f"Downloaded {len(all_observations)} observations total")

        # Convert to DataFrame
        if all_observations:
            df = pd.json_normalize(all_observations)
            if 'id' in df.columns:
                initial_count = len(df)
                # Sort so that non-null coordinates come first within each 'id'
                df = df.sort_values(by=["id", "private_location"], ascending=[True, False], na_position="last")
                # Drop duplicates, keeping the one with coordinates if available
                df = df.drop_duplicates(subset="id", keep="first")
                final_count = len(df)
                if initial_count > final_count:
                    print(f"Removed {initial_count - final_count} duplicate observations")
                print(f"Final dataset: {final_count} unique observations")
            return df
        else:
            return pd.DataFrame()

    
    def process_embedded_identifications(self, observations_df: pd.DataFrame) -> pd.DataFrame:
        """Process identifications embedded in observations data"""
        print("Processing embedded identifications...")
        
        if observations_df.empty or 'identifications' not in observations_df.columns:
            print("No identifications found in observations data")
            return pd.DataFrame()
        
        all_identifications = []
        
        for _, obs in observations_df.iterrows():
            obs_id = obs.get('id')
            identifications = obs.get('identifications', [])
            
            if isinstance(identifications, list) and identifications:
                for identification in identifications:
                    if isinstance(identification, dict):
                        # Extract identification data in format matching API response
                        processed_id = {
                            'id': identification.get('id'),
                            'uuid': identification.get('uuid'),
                            'observation.id': obs_id,
                            'user.id': identification.get('user', {}).get('id') if identification.get('user') else None,
                            'user.login': identification.get('user', {}).get('login') if identification.get('user') else None,
                            'user.name': identification.get('user', {}).get('name') if identification.get('user') else None,
                            'created_at': identification.get('created_at_details', {}).get('date'),
                            'body': identification.get('body'),
                            'category': identification.get('category'),
                            'current': identification.get('current', False),
                            'own_observation': identification.get('own_observation', False),
                            'vision': identification.get('vision', False),
                            'disagreement': identification.get('disagreement'),
                            'taxon_id': identification.get('taxon_id', 'Unknown'),
                            'hidden': identification.get('hidden', False),
                            'taxon.id': identification.get('taxon', {}).get('id') if identification.get('taxon') else None,
                            'taxon.name': identification.get('taxon', {}).get('name') if identification.get('taxon') else None,
                            'taxon.rank': identification.get('taxon', {}).get('rank') if identification.get('taxon') else None,
                            'taxon.preferred_common_name': identification.get('taxon', {}).get('preferred_common_name') if identification.get('taxon') else None,
                        }
                        all_identifications.append(processed_id)
        
        print(f"Processed {len(all_identifications)} identifications from observations")
        
        if all_identifications:
            return pd.DataFrame(all_identifications)
        else:
            return pd.DataFrame()
    
    def analyze_observations(self, observations_df: pd.DataFrame, identifications_df: pd.DataFrame, experts_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Analyze observations and categorize them into three groups:
        1. Needing review (no expert identification)
        2. Needing accuracy (expert reviewed but obscured location)
        3. Final review (expert reviewed and unobscured)
        """
        print("Analyzing observations and expert reviews...")
        
        # Create expert lookup using username (user.login)
        # Check which column contains usernames
        username_column = None
        for col in ['username', 'user_login', 'user.login']:
            if col in experts_df.columns:
                username_column = col
                break
        
        if username_column is None:
            print("Warning: No username column found in experts CSV. Expected 'username', 'user_login', or 'user.login'")
            expert_lookup = set()
        else:
            expert_lookup = set(experts_df[username_column].tolist()) if not experts_df.empty else set()
            print(f"Using '{username_column}' column for expert matching ({len(expert_lookup)} experts)")
        
        expertise_column = None
        for col in ['expertise_code', 'expertise']:
            if col in experts_df.columns:
                expertise_column = col
                break

        # Create a dictionary to track which observations have expert reviews
        expert_reviewed_obs = {}
        mismatch_experts = []  # collect mismatche expert codes
        if not identifications_df.empty:
            for _, identification in identifications_df.iterrows():
                user_login = identification.get('user.login')
                obs_id = identification.get('observation.id')
                
                if user_login in expert_lookup:
                    taxon_id_str = str(int(float(identification.get('taxon.id'))))

                    # Convert taxon_id column to string (from float or int)
                    taxon_id_series = self.taxon_name_map["taxon_id"].apply(
                        lambda x: str(int(float(x))) if pd.notna(x) else None
                    )

                    # check if identification taxon id is in the taxon name map
                    if taxon_id_str in taxon_id_series.values:
                        taxon_elcode = self.taxon_name_map.loc[taxon_id_series == taxon_id_str, "ELCODE"]
                        elcode_str = str(taxon_elcode.iloc[0])
                        user_expertise = experts_df.loc[
                            experts_df[username_column] == user_login, expertise_column
                        ]
                        user_expertise_value = str(user_expertise.iloc[0]).strip()

                        # Split and normalize expertise codes
                        if "|" in user_expertise_value:
                            expertise_codes = [code.strip().rstrip('%') for code in user_expertise_value.split('|')]
                        else:
                            expertise_codes = [user_expertise_value.rstrip('%')]

                        # Check for a match
                        match_found = any(elcode_str.startswith(code) for code in expertise_codes)

                        # confirm that the expert does qualify for expertise in this taxonomic group
                        if match_found:
                            if obs_id not in expert_reviewed_obs:
                                expert_reviewed_obs[obs_id] = []
                            
                            expert_info = experts_df.loc[
                                experts_df[username_column] == user_login
                            ].iloc[0]
                            
                            expert_reviewed_obs[obs_id].append({
                                'expert_user_login': str(user_login),
                                'expert_name': str(identification.get('user.name', 'Unknown')),
                                'expert_expertise': str(expert_info.get(expertise_column, 'Unknown')),
                                'expert_expertise_description': str(expert_info.get('expertise_description', 'Unknown')),
                                'identification_id': str(identification.get('id')),
                                'identification_taxon': str(identification.get('taxon.name', 'Unknown')),
                                'identification_date': str(identification.get('created_at')),
                                'identification_current': bool(identification.get('current', False))
                            })
                        else:
                            mismatch_experts.append({
                                'user_login': user_login,
                                'user_name': str(identification.get("user.name", "Unknown")),
                                'expertise_code': str(user_expertise.iloc[0]),
                                'taxon_elcode': elcode_str,
                                'taxon_name': identification.get('taxon.name', 'Unknown'),
                                'taxon_id': taxon_id_str
                            })
            if mismatch_experts:
                mismatch_df = pd.DataFrame(mismatch_experts)
                mismatch_df = mismatch_df.drop_duplicates(subset=["taxon_elcode", "user_login"])
                filepath = f"taxonomy/expert_mismatch_review_{self.timestamp}.csv"

                mismatch_df.to_csv(filepath, index=False)
                print(f"Saved {len(mismatch_df)} mismatched expert codes to expert_mismatch_review.csv")
            else:
                print("No mismatches found.")

        # Categorize observations
        needing_accuracy = []
        final_review = []

        for _, obs in observations_df.iterrows():
            obs_id = obs.get('id')
            taxon_id = obs.get('taxon.id')
            iNat_sciname = obs.get('taxon.name', 'Unknown')
            # Pull the taxon mapping row ONCE
            tx_row = self.taxon_name_map.loc[
                self.taxon_name_map["taxon_id"] == taxon_id
            ]
            taxon_id = obs.get('taxon.id')
            iNat_sciname = obs.get('taxon.name', 'Unknown')

            # Try 1: lookup by taxon_id
            tx_row = self.taxon_name_map.loc[
                self.taxon_name_map["taxon_id"] == taxon_id
            ]

            # Try 2: lookup by iNat name
            if tx_row.empty:
                tx_row = self.taxon_name_map.loc[
                    self.taxon_name_map["iNat_name"] == iNat_sciname
                ]
            # Try 3: fallback to binomial
            if tx_row.empty and isinstance(iNat_sciname, str):
                binomial = " ".join(iNat_sciname.split()[:2])
                tx_row = self.taxon_name_map.loc[
                    self.taxon_name_map["iNat_name"] == binomial
                ]
            
            # Try 4: check on SNAME binomial
            if tx_row.empty and isinstance(iNat_sciname, str):
                binomial = " ".join(iNat_sciname.split()[:2])
                tx_row = self.taxon_name_map.loc[
                    self.taxon_name_map["SNAME"] == binomial
                ]

            # Still nothing?
            if tx_row.empty:
                tx_values = {}
                print(f"Warning: taxon name map not found for {iNat_sciname}.")
            else:
                tx_row = tx_row.iloc[0]  
                tx_values = {
                    "SNAME": tx_row.get("SNAME"),
                    "SCOMNAME": tx_row.get("SCOMNAME"),
                    "MAJOR_GRP": tx_row.get("MAJOR_GRP"),
                    "ELCODE": tx_row.get("ELCODE"),
                    "SEID": tx_row.get("S_ELMT_ID"),
                    "TRACK": tx_row.get("TRACK"),
                    "GRANK": tx_row.get("G_RANK"),
                    "SRANK": tx_row.get("S_RANK"),
                    "RND_GRANK": tx_row.get("RNDGRANK"),
                    "FED_STATUS": tx_row.get("USESA"),
                    "OTHER_STATUS": tx_row.get("FEDSENS"),
                    "ENDEMIC": tx_row.get("ENDEMISM"),
                }
            
            # Add all observation fields plus tracking info
            base_obs_data = {
                'observation_id': obs_id,
                'observation_url': f"https://www.inaturalist.org/observations/{obs_id}",
                'observer_username': obs.get('user.login', 'Unknown'),
                'observer_fullname': obs.get('user.name', 'Unknown'),
                'observer_user_id': obs.get('user.id', 'Unknown'),
                'scientific_name': iNat_sciname,
                **tx_values,
                'description': obs.get('description', 'Unknown'),
                'observation_date': datetime.datetime.strptime(obs.get('observed_on'), '%Y-%m-%d').strftime('%Y-%m-%d'),
                'observation_date_str': f"'{str(datetime.datetime.strptime(obs.get('observed_on'), '%Y-%m-%d').strftime('%Y-%m-%d'))}",
                'created_date': obs.get('created_at'),
                'updated_date': obs.get('updated_at'),
                'quality_grade': obs.get('quality_grade', 'Unknown'),
                'location': obs.get('location', 'Unknown'),
                'positional_accuracy': obs.get('positional_accuracy'),
                'geoprivacy': obs.get('geoprivacy'),
                'obscured': obs.get('obscured', False),
                'taxon_geoprivacy': obs.get('taxon_geoprivacy'),
                'private_location': obs.get('private_location')}


            # check if the observation has been expert reviewed.
            if obs_id in expert_reviewed_obs:
                # Has expert review - check if location is obscured
                expert_reviews = expert_reviewed_obs[obs_id]

                # Add expert review information
                base_obs_data.update({
                    'expert_identification_taxon': ', '.join([er['identification_taxon'] for er in expert_reviews]),
                    'expert_identification_date': ', '.join([er['identification_date'] for er in expert_reviews]),
                    'num_expert_reviews': len(expert_reviews),
                    'all_expert_usernames': ', '.join([er['expert_user_login'] for er in expert_reviews]),
                    'all_expert_names': ', '.join([er['expert_name'] for er in expert_reviews]),
                    'expert_reviewed': True
                    
                }) 
            else:
                # No expert review
                base_obs_data.update({
                    'expert_identification_taxon': '',
                    'expert_identification_date': '',
                    'num_expert_reviews': 0,
                    'all_expert_usernames': '',
                    'all_expert_names': '',
                    'expert_reviewed': False
                })

            # Check if location is obscured
            is_obscured = (obs.get('obscured', False) or 
                            obs.get('geoprivacy') in ['obscured', 'private'] or
                            obs.get('taxon_geoprivacy') in ['obscured', 'private'])
            
            if is_obscured & pd.isna(obs.get('private_location')):
                # Mark for private coordinate enrichment
                base_obs_data['correct_coords'] = False
                needing_accuracy.append(base_obs_data)
            else:
                base_obs_data['correct_coords'] = True
                # apply corrected_latitude fields
                if not pd.isna(obs.get('private_location')):
                    base_obs_data['private_latitude'] = float(obs.get('private_location', '').split(',')[0] if obs.get('private_location') else obs.get('private_latitude', 'Not Available'))
                    base_obs_data['private_longitude'] = float(obs.get('private_location', '').split(',')[1] if obs.get('private_location') and ',' in str(obs.get('private_location')) else obs.get('private_longitude', 'Not Available'))
                    base_obs_data['corrected_latitude'] = base_obs_data['private_latitude']
                    base_obs_data['corrected_longitude'] = base_obs_data['private_longitude']
                else:
                    base_obs_data['corrected_latitude'] = float(obs.get('location', '').split(',')[0])
                    base_obs_data['corrected_longitude'] = float(obs.get('location', '').split(',')[1])
                final_review.append(base_obs_data)
        
        
        # Convert to DataFrames
        needing_accuracy_df = pd.DataFrame(needing_accuracy) if needing_accuracy else pd.DataFrame()
        final_review_df = pd.DataFrame(final_review) if final_review else pd.DataFrame()
        
        # format final_review_df for Biotics upload


        final_review_df["DIG_COM"] = final_review_df["observation_url"]

        final_review_df["UNIQUEID"] = range(1, len(final_review_df) + 1)
        final_review_df["SF_UNIQUEID"] = range(1, len(final_review_df) + 1)

        final_review_df["DISTANCE"] = final_review_df["positional_accuracy"]
        final_review_df["CONC_FEAT_"] = "Point"
        final_review_df["SF_TYPE"] = "Point"

        has_name = final_review_df["observer_fullname"].notna() & len(final_review_df["observer_fullname"])>1

        # Extract first and last names where available
        final_review_df.loc[has_name, "first_name"] = final_review_df.loc[has_name, "observer_fullname"].str.split().str[0]
        final_review_df.loc[has_name, "last_name"] = final_review_df.loc[has_name, "observer_fullname"].str.split().str[-1]

        # Fallbacks when observer_fullname is missing
        final_review_df.loc[~has_name, "first_name"] = np.nan
        final_review_df.loc[~has_name, "last_name"] = np.nan

        # --- Derived descriptive fields ---
        # DESCRIPTOR: "LastName YYYY" or "username YYYY" if name missing
        final_review_df["DESCRIPTOR"] = np.where(
            has_name,
            final_review_df["last_name"] + " " + final_review_df["observation_date"].str[:4],
            final_review_df["observer_username"] + " " + final_review_df["observation_date"].str[:4]
        )

        # LOCATOR: "iNat observation_id"
        final_review_df["LOCATOR"] = "iNat " + final_review_df["observation_id"].astype(str)

        # V_BY: "LastName, FirstName" or "observer_username" if name missing or only one part
        final_review_df["V_BY"] = np.where(
            has_name & final_review_df["observer_fullname"].str.contains(" "),
            final_review_df["last_name"] + ", " + final_review_df["first_name"],
            final_review_df["observer_username"]
        )

        # --- Location confidence fields ---
        # Convert DISTANCE to numeric (some may be strings or NaN)
        final_review_df["DISTANCE"] = pd.to_numeric(final_review_df["DISTANCE"], errors="coerce")

        # LOC_UNCERT: "Estimated" if >4.5 m, else "Negligible"
        final_review_df["LOC_UNCERT"] = np.where(
            final_review_df["DISTANCE"] > 4.5, "Estimated", "Negligible"
        )

        # REP_ACC: "High" if <50 m
        final_review_df["REP_ACC"] = np.where(
            final_review_df["DISTANCE"] > 50, "High", "Negligible"
        )

        print(f"Categorized observations:")

        print(f"  - Needing accuracy (obscured): {len(needing_accuracy_df)}")
        print(f"  - Final review (unobscured): {len(final_review_df)}")
        
        return needing_accuracy_df, final_review_df


    def run_review(self) -> None:
        """Run the complete expert review process"""
        timestamp = datetime.datetime.now().strftime("%Y%m%d")
        
        # Create output folder
        os.makedirs(self.config['output']['folder'], exist_ok=True)
        output_folder = self.config['output']['folder']
        try:
            # Load tracking list and experts
            print(f"Loaded {len(self.tracking_list)} species to track")
            experts_df = self.load_experts()
            print(f"Loaded {len(experts_df)} experts")
            
            # Download observations with both inside and outside CNHP project and zip together (required for private coords)
            observations_df = self.download_observations()
            if observations_df.empty:
                print("No observations found for tracked species. Exiting.")
                return
            
            print(f"Downloaded {len(observations_df)} observations for tracked species")
            
            # Process identifications embedded in observations data
            identifications_df = self.process_embedded_identifications(observations_df)
            
            # Analyze observations and categorize them
            needing_accuracy_df, final_review_df = self.analyze_observations(
                observations_df, identifications_df, experts_df
            )
            
            
            if not needing_accuracy_df.empty:
                needing_accuracy_file = os.path.join(output_folder, f"observations_needing_accuracy_{timestamp}.csv")
                needing_accuracy_df.to_csv(needing_accuracy_file, index=False)
                print(f"Saved {len(needing_accuracy_df)} observations needing location accuracy to: {needing_accuracy_file}")
            
            if not final_review_df.empty:
                final_review_file = os.path.join(output_folder, f"observations_for_review_{timestamp}.csv")
                final_review_df.to_csv(final_review_file, index=False)
                #self.output_shapefiles(final_review_df)
                print(f"Saved {len(final_review_df)} observations ready for review to: {final_review_file}")
                
        except Exception as e:
            print(f"Error during review process: {e}")
            import traceback
            traceback.print_exc()

def main():
    """Main entry point"""
    print("=== Offline iNaturalist Expert Review System ===")
    print(f"Started at: {datetime.datetime.now()}")
    
    reviewer = iNatScraper()
    reviewer.run_review()
    
    print(f"Completed at: {datetime.datetime.now()}")

if __name__ == "__main__":
    main()