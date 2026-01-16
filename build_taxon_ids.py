#!/usr/bin/env python3
"""
iNaturalist Taxon ID Cache Builder with Name Mapping

This script builds a dictionary mapping scientific names to taxon_ids 
for efficient batch querying of iNaturalist observations. It also tracks
how tracking list names map to iNaturalist names, handling cases where
names differ due to taxonomic updates, synonyms, or spelling differences.
Note that you may need to make manual taxon id changes, which should be logged
in a csv file in case the taxon id map needs to be rebuilt.

Usage:
    python build_taxon_cache.py

Output:
    - 'taxon_name_mappings.csv': Complete mapping table with columns:
      * iNat_name: Corresponding name found in iNaturalist (may be different)
      * taxon_id: iNaturalist taxon ID (blank if no match found)
"""

import requests
import datetime
import pandas as pd
import json
import os
import time
from typing import Dict, Optional, Tuple

class TaxonCacheBuilder:
    def __init__(self, config_file="config_private.json"):
        """Initialize with configuration file"""
        self.config_file = config_file
        self.config = self.load_config()
        
    def load_config(self) -> Dict:
        """Load configuration from JSON file"""
        
        if not os.path.exists(self.config_file):
            print(f"Config file not found: {self.config_file}")
            return {}
            
        with open(self.config_file, 'r') as f:
            return json.load(f)
    
    def get_auth_headers(self) -> Dict[str, str]:
        """Get authentication headers for API requests"""
        headers = {}
        
        # Add User-Agent (required for good API citizenship)
        if self.config.get('authentication', {}).get('user_agent'):
            headers['User-Agent'] = self.config['authentication']['user_agent']
        else:
            headers['User-Agent'] = 'iNat-TaxonCache/1.0'
        
        # Add API token if provided
        api_token = self.config.get('authentication', {}).get('api_token', '').strip()
        if api_token:
            headers['Authorization'] = f'Bearer {api_token}'
        
        return headers
    
    def preprocess_name_for_inaturalist(self, name: str) -> str:
        """
        Preprocess taxon name for iNaturalist API query by converting trinomial format
        
        Converts names like "Aster alpinus var. vierhapperi" to "Aster alpinus vierhapperi"
        which is the preferred format for iNaturalist queries.
        
        Args:
            name: Scientific name to preprocess
            
        Returns:
            Name with "var." and "ssp." removed for better iNaturalist matching
        """
        if not name or pd.isna(name):
            return name
        
        processed_name = str(name).strip()
        
        # Remove "var." and "ssp."
        processed_name = processed_name.replace(" var. ", " ")
        processed_name = processed_name.replace(" ssp. ", " ")

        
        # Clean up any double spaces
        while "  " in processed_name:
            processed_name = processed_name.replace("  ", " ")
        
        return processed_name.strip()
    
    def get_taxon_info(self, scientific_name: str) -> Tuple[Optional[int], Optional[str]]:
        """
        Get taxon_id and matched name for a scientific name using iNaturalist API
        
        Args:
            scientific_name: The scientific name to look up
            
        Returns:
            Tuple of (taxon_id, matched_name) if found, (None, None) otherwise
        """
        # Preprocess the name for better iNaturalist matching
        processed_name = self.preprocess_name_for_inaturalist(scientific_name)
        
        url = "https://api.inaturalist.org/v1/taxa"
        headers = self.get_auth_headers()
        
        params = {
            "q": processed_name,  # Use processed name for the query
            "per_page": 5,  # Get more results to find better matches
            "order": "desc",
            "order_by": "observations_count"  # Get the most observed taxon with this name
        }
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            results = data.get("results", [])
            if results:
                # Look for exact name match first
                for result in results:
                    result_name = result.get("name", "")
                    result_id = result.get("id")
                    
                    # Handle NaN values and ensure we have valid data
                    if pd.isna(result_name) or pd.isna(result_id) or not result_name or not result_id:
                        continue
                    if type(result_name)!=str:
                        continue
                    if result_name.lower() == scientific_name.lower():
                    
                        return int(result_id), result_name
                
                # If no exact match, take the first valid result (most observations)
                for result in results:
                    result_name = result.get("name", "")
                    result_id = result.get("id")
                    
                    # Handle NaN values and ensure we have valid data
                    if pd.isna(result_name) or pd.isna(result_id) or not result_name or not result_id:
                        continue
                    return int(result_id), result_name
            
            return None, None
            
        except requests.RequestException as e:
            print(f"Error looking up '{scientific_name}': {e}")
            return None, None
        except (ValueError, TypeError) as e:
            print(f"Error processing data for '{scientific_name}': {e}")
            return None, None
    
    def build_cache(self, force_rebuild: bool = False) -> Dict[str, int]:
        """
        Build taxon_id cache from tracking list and save name mappings to CSV
        
        Args:
            force_rebuild: If True, rebuild entire cache. If False, only add missing entries.
            
        Returns:
            Dictionary mapping scientific names to taxon_ids
        """
        # Load tracking list
       
        tracking_df = pd.read_csv(self.config["taxonomy"]["tracking_list"], encoding='latin1')
        
        scientific_names = tracking_df['SNAME'].unique().tolist()
        print(f"Found {len(scientific_names)} unique scientific names to process")

        # Use taxonomy folder for cache files
        taxonomy_folder = "taxonomy"
        os.makedirs("taxonomy", exist_ok=True)
        
        mapping_file = self.config["taxonomy"]["taxon_name_map"]
        out_mapping_file = os.path.join(taxonomy_folder, "taxon_name_mappings_" + datetime.datetime.now().strftime("%Y%m%d") + ".csv")
        # if the mapping file already exists, create a new one
        repeat = 1
        while os.path.exists(out_mapping_file):
            out_mapping_file = os.path.join(taxonomy_folder, "taxon_name_mappings_" + datetime.datetime.now().strftime("%Y%m%d") + "_" + str(repeat)+".csv")
            repeat +=1

        
        # Load existing mappings if they exist and we're not forcing rebuild
        if os.path.exists(mapping_file) and not force_rebuild:
            try:
                mappings_df = pd.read_csv(mapping_file)
                print(f"Loaded existing mappings with {len(mappings_df)} entries")
            except Exception as e:
                print(f"Warning: Could not load existing mappings: {e}")
        else:
            # if starting from scratch, add two columns to tracking list to store iNat_name and taxon_id
            print("rebuilding whole taxon mapping.")
            mappings_df = tracking_df
            mappings_df["taxon_id"] = None
            mappings_df["iNat_name"] = None

        matched = mappings_df[mappings_df["taxon_id"].notna()]
        print(f"Found {len(matched)} mapped taxon names.")
        # Filter out names already successfully processed (have taxon_id)
        to_match = tracking_df[~tracking_df["S_ELMT_ID"].isin(matched["S_ELMT_ID"])]
        process_number = len(to_match)
        process_status = 1
        new_rows = []

        for _, row in to_match.iterrows():
            scientific_name = row["SNAME"]
            print(f"Processing number {process_status} out of {process_number}: {scientific_name}.")

            taxon_id, matched_name = self.get_taxon_info(scientific_name)

            # Clean up taxon_id
            if taxon_id and not pd.isna(taxon_id) and str(taxon_id).lower() != "nan":
                try:
                    clean_taxon_id = int(float(taxon_id))
                except (ValueError, TypeError):
                    clean_taxon_id = None
            else:
                clean_taxon_id = None

            # Copy all existing fields from tracking_df
            new_row = row.to_dict()

            # Add / overwrite mapping fields
            new_row["taxon_id"] = clean_taxon_id
            new_row["iNat_name"] = matched_name

            new_rows.append(new_row)
            process_status += 1
            # Be kind to the API
            time.sleep(1)

        # Combine everything
        if new_rows:
            new_rows_df = pd.DataFrame(new_rows)
            mappings_df = pd.concat([mappings_df, new_rows_df], ignore_index=True)

        # Drop any duplicate SNAMEs, keeping the most recent
        mappings_df = mappings_df.drop_duplicates(subset="SNAME", keep="last")
        

        mappings_df.to_csv(out_mapping_file, index=False)
        
        return None
    
def main():
    """Main entry point"""
    print("=== iNaturalist Taxon ID Cache Builder ===")
    
    builder = TaxonCacheBuilder()
    
    print("Building taxon_id cache from tracking list...")
    print("This may take several minutes depending on the number of species.")
    print()
    
    # Build the cache (force_rebuild = True ignores any existing mapping_df and is not recommended)
    builder.build_cache(force_rebuild=False)
    
    print(f"\nDone! The cache and name mappings are ready for use with the main expert review script.")

if __name__ == "__main__":
    main()