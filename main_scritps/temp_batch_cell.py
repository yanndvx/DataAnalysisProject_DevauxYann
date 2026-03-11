# ============================================================================
# BATCH PROCESSING - ALL SUBJECTS AND SESSIONS
# ============================================================================

import os
import glob
import numpy as np

# Define base path
base_path = r'C:\Program Files\DigiMove\DigiMove\DataAnalysisProject\data\MOXY-bids'

# Get all subjects
subject_dirs = glob.glob(os.path.join(base_path, 'sub-*'))
subject_dirs = [d for d in subject_dirs if os.path.isdir(d)]

# Initialize results list
all_results = []

# Process each subject
for subject_dir in sorted(subject_dirs):
    subject_id = os.path.basename(subject_dir).split('-')[1]
    
    print(f"\n{'#'*70}")
    print(f"PROCESSING SUBJECT: sub-{subject_id}")
    print(f"{'#'*70}")
    
    # Get all sessions for this subject
    session_dirs = glob.glob(os.path.join(subject_dir, 'ses-*'))
    session_dirs = [d for d in session_dirs if os.path.isdir(d)]
    
    for session_dir in sorted(session_dirs):
        session_id = os.path.basename(session_dir).split('-')[1]
        
        print(f"\n{'-'*70}")
        print(f"Session: ses-{session_id}")
        print(f"{'-'*70}")
        
        # Check for NIRS data in this session
        nirs_dir = os.path.join(session_dir, 'nirs')
        if not os.path.exists(nirs_dir):
            print(f"  No NIRS data found for ses-{session_id}, skipping...")
            continue
        
        # Look for .snirf files to detect task type
        snirf_files = glob.glob(os.path.join(nirs_dir, f'*_recording-portalite_nirs.snirf'))
        
        if not snirf_files:
            print(f"  No .snirf file found in ses-{session_id}, skipping...")
            continue
        
        # Extract task type from filename
        snirf_filename = os.path.basename(snirf_files[0])
        if 'concentric' in snirf_filename:
            task_type = 'concentric'
            intensity_30_label = 'Con30'
            intensity_50_label = 'Con50'
            mode_name = 'Concentric'
        elif 'eccentric' in snirf_filename:
            task_type = 'eccentric'
            intensity_30_label = 'Ecc30'
            intensity_50_label = 'Ecc50'
            mode_name = 'Eccentric'
        else:
            print(f"  Unknown task type in file: {snirf_filename}, skipping...")
            continue
        
        print(f"  Task type: {task_type}")
        
        try:
            # Load NIRS data using BIDS path
            raw_path = BIDSPath(
                subject=subject_id, 
                session=session_id, 
                task=task_type,
                recording='portalite',
                suffix='nirs', 
                extension='.snirf', 
                datatype='nirs',
                root=base_path
            )
            raw = read_raw_bids(raw_path, verbose=False)
            
            # Process NIRS data
            raw_od = mne.preprocessing.nirs.optical_density(raw)
            raw_hb = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=4)
            raw_hb_filt = raw_hb.copy().filter(
                l_freq=0.01,
                h_freq=4,
                method="iir",
                iir_params=dict(order=4, ftype="butter"),
                verbose=False
            )
            
            # Get HbR channels
            ch_names = raw_hb_filt.ch_names
            hbr_channels = [ch for ch in ch_names if 'hbr' in ch.lower()]
            
            if len(hbr_channels) == 0:
                print(f"  No HbR channels found, skipping...")
                continue
            
            # Load events
            events_file = os.path.join(nirs_dir, f'sub-{subject_id}_ses-{session_id}_task-{task_type}_events.tsv')
            if not os.path.exists(events_file):
                print(f"  Events file not found: {events_file}, skipping...")
                continue
            
            events_df = pd.read_csv(events_file, sep='\t')
            
            # Process each intensity level
            window_duration = 30
            
            for intensity_label, intensity_percent in [(intensity_30_label, '30%'), (intensity_50_label, '50%')]:
                intensity_events = events_df[events_df['trial_type'] == intensity_label]
                
                if len(intensity_events) == 0:
                    print(f"    No {intensity_label} events found")
                    continue
                
                hbr_means = []
                
                # Process each series
                for idx, row in intensity_events.iterrows():
                    onset = row['onset']
                    tmin = onset
                    tmax = onset + window_duration
                    
                    try:
                        data_segment = raw_hb_filt.copy().crop(tmin=tmin, tmax=tmax)
                        data_array = data_segment.get_data()
                        
                        hbr_indices = [data_segment.ch_names.index(ch) for ch in hbr_channels if ch in data_segment.ch_names]
                        
                        if len(hbr_indices) > 0:
                            mean_hbr = np.mean(data_array[hbr_indices, :])
                            hbr_means.append(mean_hbr * 1e6)
                    except Exception as e:
                        print(f"    Error processing series at {onset}s: {e}")
                        continue
                
                if len(hbr_means) > 0:
                    mean_hbr = np.mean(hbr_means)
                    std_hbr = np.std(hbr_means)
                    
                    print(f"    {intensity_label}: {mean_hbr:.3f} ± {std_hbr:.3f} µM ({len(hbr_means)} series)")
                    
                    # Store results
                    all_results.append({
                        'Subject': subject_id,
                        'Session': session_id,
                        'Task': task_type,
                        'Contraction_Mode': mode_name,
                        'Intensity': intensity_percent,
                        'Mean_HbR_uM': mean_hbr,
                        'Std_HbR_uM': std_hbr,
                        'Number_of_Series': len(hbr_means)
                    })
                else:
                    print(f"    {intensity_label}: No valid data")
        
        except Exception as e:
            print(f"  ERROR processing session: {e}")
            continue

# Create DataFrame with all results
print(f"\n{'='*70}")
print(f"BATCH PROCESSING COMPLETE")
print(f"{'='*70}")

if len(all_results) > 0:
    results_df = pd.DataFrame(all_results)
    
    # Save to CSV
    output_path = r'C:\Program Files\DigiMove\DigiMove\DataAnalysisProject\hbr_analysis_all_subjects.csv'
    results_df.to_csv(output_path, sep=';', decimal=',', index=False)
    
    print(f"\nResults exported to: {output_path}")
    print(f"Total rows: {len(results_df)}")
    print(f"\n{results_df.to_string(index=False)}")
else:
    print("\nNo results to export.")
