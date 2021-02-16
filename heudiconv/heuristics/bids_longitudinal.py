"""
Heuristics file to be used for the Longitudinal study (data from a Stanford
GE scanner).
At the end of the classification, the function
`clean_up_unneeded_tags_from_info` simplifies the BIDS file names

Author: Pablo Velasco
Date: 2/09/2021


Notes:
Here is what I gather the different runs are:
* "FGRE", 15 slices: localizer (3 different orientations)
* "ASSET", "calibration": low-res coil calibration (is it used?)
* "efgre3d", "Inplane": inplane anatomical (anat/*_inplaneT1)
* "cni_epi": _bold
* "muxarcepi": _bold
* "'DERIVED' and 'SCREEN SAVE'": screen save
* "cni_ir_epi": _inv-<index>_IRT1 (anat)
* "cni_3dgrass", SPGR w/ varying FA: _flip-<index>_VFA (32 coils plus composite image)
"""

import os

DEFAULT_FIELDS = {
    # For the `dataset_description.json`:
    "Acknowledgements":
        "We thank Pablo Velasco (NYU Center for Brain Imaging) "
        "for preparing the BIDS dataset. "
        "TODO: whom else you want to acknowledge"
}


def create_key(subdir, file_suffix, outtype=('nii.gz', 'dicom'),
               annotation_classes=None, prefix=''):
    if not subdir:
        raise ValueError('subdir must be a valid format string')
    template = os.path.join(
        prefix,
        "{bids_subject_session_dir}",
        subdir,
        "{bids_subject_session_prefix}_%s" % file_suffix
    )
    return template, outtype, annotation_classes


def extract_task_name(prot_name):
    """ extracts task name from the protocol_name """

    prot_name = prot_name.lower()
    known_tasks = [
        'rest',
        'face',
        'conflict',
        'gamble',
        'fix',
        'films',
        'inscapes',
        'lo bold',
        'loc bold',
        'mt bold',
        'ret bold',
    ]
    for task_name in known_tasks:
        if task_name in prot_name:
            task_name = task_name.split(' bold')[0]
            break  # we keep the current value for "task"
        else:
            task_name = ''
    # ad hoc fixes:
    if task_name == 'lo':
        task_name = 'loc'

    # if we don't find a known task, try finding the keyword "task":
    if task_name == '':
        if 'task' in prot_name:
            # we want to capture what comes after "task", up to the next
            #    dash ("-") or underscore ("_"):
            task_name = prot_name.split('task')[1]
            # remove any initial "-" or "_":
            if task_name[0] == '-' or task_name[0] == '_':
                task_name = task_name[1:]
            # discard anything after the next "-" or "_":
            task_name = task_name.split('_')[0]
            task_name = task_name.split('-')[0]
            # remove any spaces we might have:
            task_name = task_name.replace(" ", "")
        else:
            task_name = 'TASK'  # fallback.  BIDS requires an alphanumeric string (no spaces...)

    return task_name


def add_series_to_info_dict(series_id, mykey, info, acq=''):
    """ adds a series to the 'info' dictionary """

    if info is None or mykey == '':
        return Exception

    if mykey in info:
        if acq == '':
            info[mykey].append({'item': series_id})
        else:
            info[mykey].append({'item': series_id, 'acq': acq})
    else:
        # if it isn't, add this key, specifying the first item:
        if acq == '':
            info[mykey] = [{'item': series_id}]
        else:
            info[mykey] = [{'item': series_id, 'acq': acq}]


def clean_up_unneeded_tags_from_info(info):
    """
    For "info" keys that only have one run, delete the "_run-<run>" tag
    """
    # Because we are going to be modifying the keys of the dictionary,
    # we get the original values now:
    orig_keys = [k for k in info.keys()]
    for k in orig_keys:
        if (
                len(info[k]) == 1
                and '_run-{item:02d}' in k[0]
        ):
            # For entries in which run number is not hard-coded:
            new_k = list(k)
            new_k[0] = new_k[0].replace('_run-{item:02d}', '')
            # remove the old k from "info" and replace it with new_k
            # while keeping the same value
            info[tuple(new_k)] = info.pop(k)
        elif '_run-01' in k[0]:
            # If "_run-01" is hardcoded, check to see if there is a
            # corresponding "_run-02". If not, remove "_run-01" from the
            # key:
            run2 = k[0].replace('_run-01', '_run-02')
            if run2 not in [kk[0] for kk in info.keys()]:
                new_k = list(k)
                new_k[0] = new_k[0].replace('_run-01', '')
                # remove the old k from "info" and replace it with new_k
                # while keeping the same value
                info[tuple(new_k)] = info.pop(k)


def sort_seqinfo_by_series_date_time(seqinfo):
    """
    Sorts the 'seqinfo' list according to the series acquisition date
    and time (using the "series_uid" to make sure there are no repeats)

    This assures that all the entries are in acquisition order (which
    is not always the case).

    input: seqinfo: original seqinfo list
    output: sortedSeqinfo: seqinfo sorted list
    """
    # sort by concatenated date and time (strings):
    # (use set to only keep unique ones):
    dateTimes = set([s.date + s.time for s in seqinfo if s.date and s.time])
    sortedSeqinfo = []
    for dt in sorted(dateTimes):
        dTseries = [
            ss for ss in seqinfo if (
                    ss.date
                    and ss.time
                    and ss.date + ss.time == dt
            )
        ]
        # sort series with identical date and time by series_uid:
        for sid in sorted([s.series_uid for s in dTseries if s.series_uid]):
            for ss in dTseries:
                if ss.series_uid == sid:
                    sortedSeqinfo.append(ss)

        # Now, add the series which do not have series_uid:
        for ss in dTseries:
            if ss not in sortedSeqinfo:
                sortedSeqinfo.append(ss)

    # Now, add the series which do not have date or time:
    for ss in seqinfo:
        if ss not in sortedSeqinfo:
            sortedSeqinfo.append(ss)

    return sortedSeqinfo


def infotodict(seqinfo):
    """
    Heuristic evaluator for determining which runs belong where allowed
    template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    ###   NIFTI and DICOM   ###
    # anat:
    t1 = create_key('anat', 'acq-{acq}_run-{item:02d}_T1w')
    t2 = create_key('anat', 'acq-{acq}_run-{item:02d}_T2w')
    inplaneT1 = create_key('anat', 'run-{item:02d}_inplaneT1')
    pd_bias_body = create_key('anat', 'acq-biasBody_run-{item:02d}_PD')
    pd_bias_receive = create_key('anat', 'acq-biasReceive_run-{item:02d}_PD')
    ir = create_key('anat', 'inv-{item:02d}_IRT1')
    variable_fa = create_key('anat', 'flip-{item:02d}_VFA')

    # func:
    # For the functional runs, we want to have a different key for each
    # task. Since we don't know a-priori what the task names will be, we
    # don't create the different keys here, but we'll create them
    # dynamically, as needed.

    # diffusion:
    dwi = create_key('dwi', 'acq-{acq}_run-{item:02d}_dwi')
    dwi_sbref = create_key('dwi', 'acq-{acq}_run-{item:02d}_sbref')

    # fmaps:
    # For the SE-EPI field-map runs (for topup), both for fmri and dwi,
    # we want to have a different key for each PE direction.  Rather
    # than creating the different keys here, we'll create them
    # dynamically, as needed.
    fmap_gre_mag = create_key('fmap', 'acq-GRE_run-{item:02d}_magnitude')  # GRE fmap
    fmap_gre_phase = create_key('fmap', 'acq-GRE_run-{item:02d}_phasediff')  #

    ###  DICOM only   ###
    # These are images we don't want for analysis, but we still want to
    # keep a copy of the DICOMs in the 'sourcedata' folder.  We manage
    # that by specifying the 'outtype' to be only 'dicom':
    # anat:
    t2_scout = create_key('anat', 'acq-Scout_run-{item:02d}_T2w', outtype=('dicom',))
    # Misc:
    screen_save = create_key('misc', 'run-{item:02d}_screensave', outtype=('dicom',))

    info = {t1: [], t2: [], inplaneT1: [], pd_bias_body: [], pd_bias_receive: [],
            dwi: [], dwi_sbref: [],
            fmap_gre_mag: [], fmap_gre_phase: [],
            t2_scout: [], screen_save: []}
    run_numbers = {}  # dict with the run number for each task

    # because not always are series sorted by acquisition date/time,
    # sort them:
    seqinfo = sort_seqinfo_by_series_date_time(seqinfo)

    for idx, s in enumerate(seqinfo):
        # s is a namedtuple with fields equal to the names of the columns
        # found in the dicominfo.tsv file

        # we don't need the original case of the protocol name:
        prot_name = s.protocol_name.lower()
        acq = ''

        ###   T1w   ###
        # 1) High resolution T1w:
        # single volume, protocol name including T1, MPRAGE, MP-RAGE, MPR,...
        if (
            s.dim4 == 1
            and 'fl' in s.sequence_name
            and (
                't1' in prot_name
                or ('mp' in prot_name and 'rage' in prot_name)
                or 'mpr' in prot_name
            )
        ):
            acq = 'highres' + direction  # if direction is empty, aqc='highres'

            # If this image is NOT normalized, check if another
            # series has identical acquisition date and time.
            # If so, we'll keep only the normalized version, and for
            # this one we keep just the DICOM.
            # (Note: older versions of heudiconv don't include 'time'):
            # Otherwise, we extract it:
            if (
                'NORM' not in s.image_type
                and another_series_has_identical_start_time(seqinfo, idx)
            ):
                info[t1_dicom].append({'item': s.series_id, 'acq': acq})
            else:
                info[t1].append({'item': s.series_id, 'acq': acq})

        # 2) High-res inplane:
        # single volume, 'efgre3d' sequence, 'inplane'
        elif (
            s.dim4 == 1
            and 'efgre3d' in s.sequence_name
            and 'inplane' in s.series_description.lower()
        ):
            info[inplaneT1].append({'item': s.series_id})

        # 3) FSE T1w:
        # single volume, series description includes TSE or FSE,
        # protocol name includes T1
        elif (
            s.dim4 == 1
            and 't1' in prot_name
            and 'tse' in s.sequence_name
        ):
            acq = 'fse' + direction  # if direction is empty, aqc='fse'
            info[t1].append({'item': s.series_id, 'acq': acq})

        ###   T2w   ###
        # 1) Scout:
        if s.sequence_name == 'fgre':
            info[t2_scout].append({'item': s.series_id})

        # 2) Standard high-res T2w used for cortical segmentation:
        # single volume, protocol name including T2, T2w, TSE, SPACE, SPC:
        elif (
                s.dim4 == 1
                and (
                        'T2' in s.protocol_name
                        or 'tse' in prot_name
                        or 'space' in prot_name
                        or 'spc' in prot_name
                )
        ):
            acq = 'highres' + direction  # note: if direction is empty, aqc='highres'

            # If this image is NOT normalized, check if another
            # series has identical acquisition date and time.
            # If so, we'll keep only the normalized version, and for
            # this one we keep just the DICOM.
            # (Note: older versions of heudiconv don't include 'time'):
            # Otherwise, we extract it:
            if (
                    'NORM' not in s.image_type
                    and another_series_has_identical_start_time(seqinfo, idx)
            ):
                info[t2_dicom].append({'item': s.series_id, 'acq': acq})
            else:
                info[t2].append({'item': s.series_id, 'acq': acq})

        # 3) Fast Spin-Echo used as a localizer:
        # single volume, sequence name: 'h2d1' ('haste')"
        elif (
                s.dim4 == 1
                and 'h2d1' in s.sequence_name
        ):
            info[t2_dicom].append({'item': s.series_id, 'acq': 'haste'})

        ###   PD   ###
        # BIAS images  (for coil sensitivity estimation) are typically
        # PD-weighted
        elif (
                'ASSET calibration' in s.series_description
                and 'efgre3d' in s.sequence_name
        ):
            info[pd_bias_receive].append({'item': s.series_id})

        ###   INVERSION_RECOVERY   ###
        elif 'cni_ir_epi' in s.sequence_name:
            # We could parse the TI from the series_description, but
            # we don't need to. It will be in the json file
            add_series_to_info_dict(s.series_id, ir, info)

        ###   VARIABLE FLIP ANGLE   ###
        elif (
                'cni_3dgrass' in s.sequence_name
                and 'deg' in s.series_description
        ):
            # We could parse the FA from the series_description, but
            # we don't need to. It will be in the json file
            add_series_to_info_dict(s.series_id, variable_fa, info)

        ###   FUNCTIONAL   ###
        elif (
                s.dim4 >= 4
                and (
                        'cni_epi' in s.sequence_name
                        or 'muxarcepi' in s.sequence_name
                )
                and not 'DERIVED' in s.image_type
        ):

            # check task name:
            task = extract_task_name(s.series_description)
            # find the run number for this task:
            if run_numbers.get(task):
                run_numbers[task] += 1
            else:
                run_numbers[task] = 1

            # dictionary key specific for this task type:
            mykey = create_key(
                'func',
                'task-%s_run-%02d_bold' % (task, run_numbers[task])
            )
            add_series_to_info_dict(s.series_id, mykey, info)

        ###   FIELD MAPS   ###

        # A) Spin-Echo distortion maps to be used with topup:
        #  (note: we only look for magnitude images, because we'll catch
        #   the phase images below
        #   Also, we exclude _sbref because the sbref from MB diffusion
        #   also have epse2d in the sequence_name
        elif (
                s.dim4 <= 3
                and 'epse2d' in s.sequence_name
                and 'M' in s.image_type
                and (
                        'dist' in prot_name
                        or 'map' in prot_name
                        or 'field' in prot_name
                )
                and not s.series_description.lower().endswith('_sbref')
        ):
            # check PE direction:
            run_dir = direction or 'normal'

            # is phase image present?
            # At least for Siemens systems, if magnitude/phase was
            # selected, the phase images come as a separate series
            # immediatelly following the magnitude series.
            # (note: make sure you don't check beyond the number of
            # elements in seqinfo...)
            if (
                    idx + 1 < len(seqinfo)
                    and seqinfo[idx + 1].protocol_name == s.protocol_name
                    and 'P' in seqinfo[idx + 1].image_type
            ):
                # we have a magnitude/phase pair:

                # dictionary keys specific for this SE-fmap direction:
                mykey_mag = create_key(
                    'fmap',
                    'acq-fMRI_rec-magnitude_dir-%s_run-{item:02d}_epi' % run_dir
                )
                mykey_pha = create_key(
                    'fmap',
                    'acq-fMRI_rec-phase_dir-%s_run-{item:02d}_epi' % run_dir
                )
                add_series_to_info_dict(s.series_id, mykey_mag, info)
                add_series_to_info_dict(
                    seqinfo[idx + 1].series_id, mykey_pha, info
                )

            else:
                # we only have a magnitude image

                # dictionary key specific for this SE-fmap direction:
                mykey = create_key(
                    'fmap',
                    'acq-fMRI_dir-%s_run-{item:02d}_epi' % run_dir
                )
                add_series_to_info_dict(s.series_id, mykey, info)

        # B) GRE fmap:
        elif (
                'fm2d' in s.sequence_name
                and ('fieldmap' in prot_name or 'field_map' in prot_name)
        ):
            if s.image_type[2] == 'M':
                # magnitude image
                info[fmap_gre_mag].append({'item': s.series_id})
            elif s.image_type[2] == 'P':
                # phase image:
                info[fmap_gre_phase].append({'item': s.series_id})


        ###   DIFFUSION   ###

        # We could also check: (s.image_type[2] == 'DIFFUSION')
        elif 'ep_b' in s.sequence_name:
            # Siemens product diffusion sequence

            # This is not very rigorous, but in general, diffusion runs
            # will have more than a couple of volumes.  I'll use 5, to
            # be safe:
            if s.dim4 >= 5:
                # this is a standard diffusion acquisition
                acq = str(s.dim4) + 'vols'
                info[dwi].append({'item': s.series_id, 'acq': acq})

                # check to see if the previous run is a SBREF:
                if (
                        previous_series_is_sbref(seqinfo, idx)
                        and 'epse2d' in seqinfo[idx - 1].sequence_name
                ):
                    info[dwi_sbref].append(
                        {'item': seqinfo[idx - 1].series_id, 'acq': acq}
                    )

            else:
                # this is a fmap for diffusion.

                # Get the PE direction, for topup:
                run_dir = direction or 'normal'

                # dictionary key specific for this SE-fmap direction:
                mykey = create_key(
                    'fmap',
                    'acq-dwi_dir-%s_run-{item:02d}_epi' % run_dir
                )
                add_series_to_info_dict(s.series_id, mykey, info)

                if previous_series_is_sbref(seqinfo, idx):
                    # TO-DO: for now, extract the _sbref dwi fmap image
                    # as DICOMs, because BIDS doesn't allow them.
                    mykey = create_key(
                        'fmap',
                        'acq-dwi_dir-%s_run-{item:02d}_sbref' % run_dir,
                        outtype=('dicom',)
                    )
                    add_series_to_info_dict(
                        seqinfo[idx - 1].series_id, mykey, info
                    )

        ###   SCREEN SAVE   ###

        elif (
            'DERIVED' in s.image_type
            and 'SCREEN SAVE' in s.image_type
        ):
            info[screen_save].append({'item': s.series_id})

    clean_up_unneeded_tags_from_info(info)
    return info
