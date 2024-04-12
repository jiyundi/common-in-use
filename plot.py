import matplotlib.pyplot as plt

fig = plt.figure(figsize=(7,9),dpi=100) # height=7
plt.subplots_adjust(hspace=0.15, wspace=0.20) # h=height
gs = fig.add_gridspec(2, 1,
                      height_ratios=[1,1],
                      width_ratios=[1])
ax_orig_spec = fig.add_subplot(gs[0, :])
ax_orig_erro = fig.add_subplot(gs[1, :])

ax_orig_spec.plot(arr_obs[0],
                  arr_obs[1], 
                  label='Observed spectrum',
                  color='black', linewidth=1, zorder=11)
ax_orig_erro.plot(arr_obs[0],
                  arr_obs[2], 
                  color='black', linewidth=1, zorder=11)

ax_orig_spec.set_ylabel('Flux')
ax_orig_erro.set_ylabel('Error in Flux')
ax_orig_spec.minorticks_on()
ax_orig_erro.minorticks_on()
ax_orig_spec.grid(alpha=0.75, linestyle='--', zorder=-1)
ax_orig_erro.grid(alpha=0.75, linestyle='--', zorder=-1)
ax_orig_spec.legend(loc='upper left', prop={'size': 8})
                      
plt.savefig(spec1dfilepath+".jpg", dpi=150, bbox_inches='tight')
