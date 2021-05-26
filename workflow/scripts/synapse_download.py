import synapseclient
import synapseutils


token = open(snakemake.input['token']).readline().strip()
syn = synapseclient.Synapse()
syn.login(authToken=token)

files = synapseutils.syncFromSynapse(syn, snakemake.params['project_id'],
                                     path=snakemake.output['lrgasp_dir'])
