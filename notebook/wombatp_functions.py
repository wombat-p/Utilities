    
def DataDir(wdgt):
    data_dir = wdgt.Text(
        placeholder='Type something',
        value = './data_small',
        description='',
        disabled=False,
    )
    
    mode = wdgt.RadioButtons(
    options=['Proline', 'Transproteomic', 'Compomic', "OpenMS"],
    value='Proline',
    continuous_update=True,
    disabled=False
)

    print("\nPath to your data directory.")
    print('This directory will be used to infer further files locations.')
    display(data_dir)
    
    print('Select the mod')
    display(mode)
    return [data_dir,mode]



def Arguments(wdgt, pd, data_dir, selection):
    path = 'unimod_searchgui_mapping.txt'
    df = pd.read_csv("unimod_searchgui_mapping.txt", sep='\t')
    modification = df["Name"].values.tolist()
    dictionary = {}
    
    run_statistics = wdgt.Checkbox(
            value=False,
            description='',
            disabled=False,
            indent=False,
            continuous_update=True,
        )
    
    dictionary["run_statistics"] = run_statistics
    if "Proline" == selection.value or "Compomic" == selection.value or "Transproteomic" == selection.value:
        raws = wdgt.Text(
            placeholder='Type something',
            description='',
            disabled=False,
            value=data_dir.value + '/*.raw'
        )
        dictionary["raws"] = raws

        fasta = wdgt.Text(
            placeholder='Type something',
            description='',
            disabled=False,
            value=data_dir.value + '/*.fasta'
        )
        dictionary["fasta"] = fasta
        
        precursor_mass_tolerance = wdgt.FloatSlider(
            value=10,
            min=0,
            max=100.0,
            step=0.1,
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='.1f',
        )
        dictionary["precursor_mass_tolerance"] = precursor_mass_tolerance   
        
        fragment_mass_tolerance = wdgt.FloatSlider(
            value=0.5,
            min=0,
            max=2.0,
            step=0.01,
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='.2f',
        )
        dictionary["fragment_mass_tolerance"] = fragment_mass_tolerance
    
        miscleavages = wdgt.IntSlider(
            value=2,
            min=0,
            max=10,
            step=1,
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='d'
        )
        dictionary["miscleavages"] = miscleavages
    
        AcetylationIndex = modification.index("Acetylation of protein N-term")
        OxidationIndex = modification.index("Oxidation of M")
        variable_mods = wdgt.SelectMultiple(
            options=modification,
            value=(modification[AcetylationIndex] , modification[OxidationIndex]), # FIXME: should use names instead
            rows=3,
            disabled=False,
            continuous_update=True,
            layout = wdgt.Layout(height="250px",width="300px")
        )
        dictionary["variable_mods"] = variable_mods
             
        experiment_design = wdgt.Text(
            placeholder='Type something',
            description='',
            disabled=False,
            value=data_dir.value + '/pxd001819.txt'
        )
        dictionary["experiment_design"] = experiment_design
        
        if "Proline" == selection.value:
            lfq_param = wdgt.Text(
                placeholder='Type something',
                description='',
                disabled=False,
                value=data_dir.value + '/lfq_param_file.txt'
            )
            dictionary["lfq_param"] = lfq_param
        
        if "Compomic" == selection.value or "Transproteomic" == selection.value :
            enzyme = wdgt.Text(
                placeholder='Type something',
                description='',
                disabled=False,
                value = "Trypsin (no P rule)"
            )
            dictionary["enzyme"] = enzyme
        
        if "Transproteomic" == selection.value :
            fdr_peptide_threshold = wdgt.FloatSlider(
                value=0.5,
                min=0,
                max=2.0,
                step=0.01,
                disabled=False,
                continuous_update=True,
                orientation='horizontal',
                readout=True,
                readout_format='.2f',
            )
            dictionary["fdr_peptide_threshold"] = fdr_peptide_threshold
    
            quantification_fdr = wdgt.FloatSlider(
                value=0.5,
                min=0,
                max=2.0,
                step=0.01,
                disabled=False,
                continuous_update=True,
                orientation='horizontal',
                readout=True,
                readout_format='.2f',
            )
            dictionary["quantification_fdr"] = quantification_fdr
        
        print("\nPath to input data.(must be surrounded with quotes)")
        print('If you have several paths, fill in the source folder like: "RAWFOLDER/*.raw" \nelse: "RAWFOLDER/FileName.raw"')
        display(raws)
    
        print("\nFasta file for database search")
        print('If you have several paths, fill in the source folder like: "../data/*.fasta" \nelse: "../data/yeast_UPS.fasta"')
        display(fasta)
        
        if "Proline" == selection.value:
            print("\nParameter file for Proline")
            display(lfq_param)

        print("\nMass tolerance of precursor mass (ppm)")
        display(precursor_mass_tolerance)
    
        print("\nMass tolerance of fragment mass bin (Da)")
        display(fragment_mass_tolerance)
    
        print("\nNumber of allowed miscleavages")
        display(miscleavages)

        print("\nVariable modifications ('Oxidation of M', see Search modifications)")
        display(variable_mods) 
        
        if "Compomic" == selection.value or "Transproteomic" == selection.value :
            print("\nEnzymatic cleavage (e.g. 'Trypsin', see SearchGUI enzymes)")
            display(enzyme)

        print("text-file containing 2 columns: first with mzDB file names and second with names for experimental conditions")
        display(experiment_design)
       
    if "Transproteomic" == selection.value :
        print("\nFalse Discovery Rate Peptide threshold")
        display(fdr_peptide_threshold)
        print("\nFalse Discovery Rate quantification")
        display(quantification_fdr)
    
    print("\nRun statistics :")
    display(run_statistics)
    return dictionary




def CleanDir(wdgt, os):
    def clr(self):
        os.system("bash clean_dir.sh")
    clean = wdgt.Button(
        description='Clean directory',
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon=''
    )
    display(clean)
    clean.on_click(clr)
    

    
#def GenerateCmdLine(wdgt, workflows_dir, multiprocessing, psutil, output, selection): #TODO
def GenerateCmdLine(wdgt, multiprocessing, psutil, output, selection):
    button = wdgt.Button(
        description='Generate',
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon=''
    )

    cmd = wdgt.Textarea(
        value=selection.value,
        placeholder='',
        description='',
        disabled=True,
        layout = wdgt.Layout(height="75px",width="auto")
    )
    if selection.value != "OpenMS":
        variable_mods = ",".join(output["variable_mods"].value)
    def GenerateCommandLine(cdm):
        memory = int((psutil.virtual_memory().available/1000000000)-5)
        cpu = multiprocessing.cpu_count()
        if selection.value == "Compomic":
            cmd.value = 'Compomics-Workflow/Nextflow/./nextflow run Compomics-Workflow/Nextflow/main.nf --raws "{0}" --fasta "{1}" --miscleavages {2} \
--fragment_mass_tolerance {3} --precursor_mass_tolerance {4} --enzyme "{5}" \
--variable_mods "{6}" --experiment_design "{7}" --run_statistics {8} --max_cpus {9} --max_memory {10}GB \
-profile docker -with-report -with-trace -with-timeline'.format(    output["raws"].value,
                                                                    output["fasta"].value,
                                                                    output["miscleavages"].value,
                                                                    output["fragment_mass_tolerance"].value,
                                                                    output["precursor_mass_tolerance"].value,
                                                                    output["enzyme"].value,
                                                                    variable_mods,
                                                                    output["experiment_design"].value,
                                                                    output["run_statistics"].value,
                                                                    cpu,
                                                                    memory
                                                               )                                                                
            
        if selection.value == "Proline":
            cmd.value = 'Proline-Workflow/Nextflow/./nextflow run Proline-Workflow/Nextflow/main.nf --raws "{0}" --fasta "{1}" --precursor_mass_tolerance {2} \
--fragment_mass_tolerance {3} --miscleavages {4} --variable_mods "{5}" \
--experiment_design "{6}" --lfq_param "{7}" \
 --run_statistics {8} --max_cpus {9} --max_memory {10}GB -profile docker -with-report -with-trace -with-timeline'.format(     output["raws"].value,
                                                                                                       output["fasta"].value, 
                                                                                                       output["precursor_mass_tolerance"].value, 
                                                                                                       output["fragment_mass_tolerance"].value,
                                                                                                       output["miscleavages"].value,
                                                                                                       variable_mods, 
                                                                                                       output["experiment_design"].value,
                                                                                                       output["lfq_param"].value,
                                                                                                       output["run_statistics"].value,
                                                                                                       cpu,
                                                                                                       memory)
        if selection.value == "Transproteomic":
            cmd.value ='./Transproteomic-Pipeline/Nextflow/nextflow run ./Transproteomic-Pipeline/Nextflow/main.nf --raws "{0}" --fasta "{1}" --miscleavages {2} --fragment_mass_tolerance {3} \
--precursor_mass_tolerance {4} --enzyme "{5}" --variable_mods {6} \
--fdr_peptide_threshold {7} --quantification_fdr {8} --experiment_design {9} --max_cpus {10} --max_memory \
{11}GB -profile docker -with-report -with-trace -with-timeline'.format(     output["raws"].value,
                                                                            output["fasta"].value,
                                                                            output["miscleavages"].value,
                                                                            output["fragment_mass_tolerance"].value,
                                                                            output["precursor_mass_tolerance"].value,
                                                                            output["enzyme"].value,
                                                                            variable_mods,
                                                                            output["fdr_peptide_threshold"].value,
                                                                            output["quantification_fdr"].value,
                                                                            output["experiment_design"].value,
                                                                            cpu,
                                                                            memory
                                                                   )
            
        if selection.value == "OpenMS":
            cmd.value = "./OpenMS-ProteomicsLFQ/Nextflow/./nextflow run nf-core/proteomicslfq -r dev -profile docker --input https://raw.githubusercontent.com/bigbio/proteomics-metadata-standard/master/annotated-projects/PXD001819/PXD001819.sdrf.tsv \
--database https://raw.githubusercontent.com/wombat-p/Transproteomic-Pipeline/dev/Results/yeast_UPS.fasta --add_decoys -with-report -with-trace -with-timeline"


    
    print("Generated command line :")
    display(button)
    button.on_click(GenerateCommandLine)
    return cmd

def Run(wdgt, os, output, mode):
    run = wdgt.Button(
        description='Run',
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon=''
    )
    
    def RunButton(self):
        os.system(output["cmd"].value)
    
    display(run)
    run.on_click(RunButton)

def Resume(wdgt, os, output):
    resume = wdgt.Button(
        description='Resume',
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon=''
    )
    def resumeButton(self):
         os.system(output["cmd"].value + ' -resume')

    display(resume)
    resume.on_click(resumeButton)