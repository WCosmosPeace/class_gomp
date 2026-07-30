#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import yaml
from getdist import loadMCSamples, plots


GOMP_COLOR = '#4477AA'
DM_DECAY_COLOR = '#CC6677'
LINE_WIDTH = 1.5
CHAIN_USAGE = (
    'Usage: python tria.py <DM-decay-chain>/mcmc <Gomp-chain>/mcmc\n'
    'Example: python tria.py Gomp_mnu_decay/mcmc Gomp_mnu/mcmc'
)

PLOT_PARAMETERS = [
    # Core cosmology
    'n_s',
    'H0',
    'theta_s_100',
    'Omega_m',
    'Omega_b',
    'omega_cdm',
    # Primordial amplitude
    'A_s',
    'S8',
    # Reionization
    'tau_reio',
    # Neutrino mass
    'mnu',
    # Dark-matter and Gompertz models
    'decay_rate',
    'alpha_gomp',
    'beta_gomp',
    # Nuisance parameters
    'A_act',
    'P_act',
    'A_planck',
    'YHe',
    # Geometry
    'rs_drag',
    'age',
    # Fit diagnostics
    'chi2',
    'chi2__lowlEE',
    'chi2__lowlTT',
    'chi2__Planck',
    'chi2__ACT',
    'chi2__BAO',
    'chi2__QDW',
]

# ACT DR6 Table 5, P-ACT-LB column.
REFERENCE_VALUES = {
    'A_s': 2.1327557162e-9,  # 1e-10 * exp(3.060)
    'n_s': 0.9743,
    'theta_s_100': 1.04086,
    'omega_cdm': 0.1179,
    'Omega_m': 0.3032,
    'Omega_b': 0.04847,
    'H0': 68.22,
    'S8': 0.8169,
    'tau_reio': 0.0632,
    'YHe': 0.245937,
    'age': 13.772,
    'rs_drag': 147.45,
}


def apply_plot_style():
    plt.rcParams.update({
        'lines.dash_joinstyle': 'round',
        'lines.dash_capstyle': 'round',
        'lines.solid_joinstyle': 'round',
        'lines.solid_capstyle': 'round',
        'font.family': 'serif',
        'font.serif': ['Palatino'],
        'text.usetex': True,
        'xtick.top': True,
        'xtick.direction': 'in',
        'ytick.right': True,
        'ytick.direction': 'in',
        'legend.frameon': False,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        'savefig.transparent': True,
    })


def prepare_sample(chain_root):
    sample = loadMCSamples(chain_root, settings={'ignore_rows': 0.2})
    parameter_names = sample.getParamNames()

    # Normalize names and labels so both chain sets use the same plot axes.
    for parameter in parameter_names.names:
        if parameter.name == 'm_ncdm':
            parameter.name = 'mnu'
            parameter.label = r'\sum m_\nu'
        elif parameter.name == 'alpha_gomp':
            # alpha_gomp already stores ln(alpha) in these chains.
            parameter.label = r'\ln \alpha'

    sample.setParamNames(parameter_names)
    sample.updateBaseStatistics()

    # Add S8 only when it is not already stored in the chain.
    if 'S8' not in parameter_names.list():
        try:
            params = sample.getParams()
            sample.addDerived(
                params.sigma8 * (params.Omega_m / 0.3) ** 0.5,
                name='S8',
                label='S_8',
            )
            sample.updateBaseStatistics()
            print(f'>>> Added derived parameter S8 for {chain_root}')
        except Exception as error:
            print(f'>>> Could not add S8 for {chain_root}: {error}')

    return sample


def main(chain_roots):
    if len(chain_roots) != 2:
        raise SystemExit(f'Expected exactly two chains.\n{CHAIN_USAGE}')

    apply_plot_style()

    samples = []
    labels = []
    line_args = []
    dm_mass_ev = None

    for index, chain_root in enumerate(chain_roots):
        sample = prepare_sample(chain_root)
        has_dm_decay = 'decay_rate' in sample.getParamNames().list()

        # The first chain defines the plotted parameters and title statistics.
        # It must contain decay_rate; the second chain must be the Gomp model.
        expected_dm_decay = index == 0
        if has_dm_decay != expected_dm_decay:
            raise SystemExit(f'Incorrect chain order.\n{CHAIN_USAGE}')

        if has_dm_decay:
            # Read the fixed DM mass from the configuration used for this chain.
            input_path = Path(f'{chain_root}.input.yaml')
            with input_path.open(encoding='utf-8') as stream:
                config = yaml.safe_load(stream) or {}

            mass_config = config.get('params', {}).get('mDM')
            if not isinstance(mass_config, dict) or 'value' not in mass_config:
                raise ValueError(
                    f'Could not find params.mDM.value in {input_path}'
                )

            chain_mass_ev = float(mass_config['value'])
            if chain_mass_ev <= 0:
                raise ValueError(
                    f'params.mDM.value must be positive in {input_path}: '
                    f'{chain_mass_ev}'
                )
            if dm_mass_ev is not None and chain_mass_ev != dm_mass_ev:
                raise ValueError(
                    'DM-decay chains use inconsistent mDM values: '
                    f'{dm_mass_ev:g} eV and {chain_mass_ev:g} eV'
                )

            dm_mass_ev = chain_mass_ev
            label = r'Gomp $+\,\sum m_\nu+$ DM decay'
            # Draw the red DM-decay contours above the blue Gomp contours.
            style = {
                'color': DM_DECAY_COLOR,
                'ls': '--',
                'lw': LINE_WIDTH,
                'zorder': 3,
            }
            print(f'>>> Read mDM = {dm_mass_ev:g} eV from {input_path}')
        else:
            label = r'Gomp $+\,\sum m_\nu$'
            style = {
                'color': GOMP_COLOR,
                'ls': '-',
                'lw': LINE_WIDTH,
                'zorder': 2,
            }

        samples.append(sample)
        labels.append(label)
        line_args.append(style)

    # Samples stay in [DM decay, Gomp] order so decay_rate is available,
    # while the legend is displayed in the more natural [Gomp, DM decay] order.
    legend_order = [1, 0]

    # The first chain defines the axes and the title statistics. 
    # Put the DM-decay chain first when decay_rate must be included.
    available = set(samples[0].getParamNames().list())
    parameters = [
        name for name in PLOT_PARAMETERS if name in available
    ]
    missing = [
        name for name in PLOT_PARAMETERS if name not in available
    ]

    print('\n=== Parameters to plot ===')
    print(parameters)
    if missing:
        print('\n=== Missing from the reference chain; skipped ===')
        print(missing)

    plotter = plots.get_subplot_plotter(subplot_size=0.9)
    plotter.settings.linewidth = LINE_WIDTH
    plotter.settings.axes_fontsize = 9
    plotter.settings.lab_fontsize = 11
    plotter.settings.legend_fontsize = 18
    plotter.settings.title_limit = 1    # Show the first chain's 68% interval above each diagonal panel.
    plotter.settings.title_limit_fontsize = 10

    plotter.triangle_plot(
        samples,
        parameters,
        filled=False,
        legend_labels=labels,
        contour_colors=[style['color'] for style in line_args],
        contour_ls=[style['ls'] for style in line_args],
        contour_lws=[style['lw'] for style in line_args],
        contour_args=[
            {'zorder': style['zorder']} for style in line_args
        ],
        line_args=line_args,
        label_order=legend_order,
        legend_ncol=1,
        legend_loc='upper right',
        # Dashed reference lines use the P-ACT-LB values listed above.
        markers=REFERENCE_VALUES,
    )

    # Match the reordered legend text to the blue Gomp and red DM-decay lines.
    for text, color in zip(
        plotter.legend.get_texts(),
        (GOMP_COLOR, DM_DECAY_COLOR),
    ):
        text.set_color(color)

    if dm_mass_ev is not None:
        # Use the largest unit that keeps the displayed mass at least one.
        mass_value = dm_mass_ev
        mass_unit = 'eV'
        for unit, scale in (
            ('TeV', 1e12),
            ('GeV', 1e9),
            ('MeV', 1e6),
            ('keV', 1e3),
        ):
            if dm_mass_ev >= scale:
                mass_value = dm_mass_ev / scale
                mass_unit = unit
                break

        # A canvas draw is required before the legend size and position exist.
        plotter.fig.canvas.draw()
        renderer = plotter.fig.canvas.get_renderer()
        legend_box = plotter.legend.get_window_extent(renderer)
        legend_box = legend_box.transformed(
            plotter.fig.transFigure.inverted()
        )
        legend_fontsize = plotter.legend.get_texts()[0].get_fontsize()
        plotter.fig.text(
            legend_box.x1,
            legend_box.y0 - 0.003,
            rf'$m_\chi = {mass_value:g}\,\mathrm{{{mass_unit}}}$',
            color=DM_DECAY_COLOR,
            fontsize=legend_fontsize,
            ha='right',
            va='top',
        )

        mass_tag = (
            f'{mass_value:g}'.replace('.', 'p').replace('+', '')
            + mass_unit
        )
        output_path = Path(f'triangle_{mass_tag}.pdf')
    else:
        output_path = Path('triangle_Gomp_mnu_comparison.pdf')

    plotter.export(str(output_path))
    print(f'\nSaved: {output_path.resolve()}')


if __name__ == '__main__':
    main(sys.argv[1:])
