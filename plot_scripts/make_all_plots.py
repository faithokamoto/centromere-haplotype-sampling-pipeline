## == File set-up == ##
# ARIAL

import argparse
import random

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np

# Basic plotting settings
plt.style.use('BME163')
NUM_STEPS = 101
POINT = {'marker': 'o', 'markeredgewidth': 0, 'linewidth': 0, 'markersize': 0.5}
POINT_OVERLAP_DIST = POINT['markersize'] / 75

## == Plotting parameters == ##

# Total figure
figure_size = (6, 2.5) # (width, height)
figure_dpi = 600
# Panels
main_panel_loc = (0.5, 0.375, 4.5, 1.5) # (left, bottom, width, height)
right_panel_loc = (5.1, 0.375, 0.2, 1.5)
main_panel_lim = {'x': (0, 8), 'y': (75, 100)}
right_panel_lim = (0, NUM_STEPS)
right_panel_ticks = (8, 9, 10, 11, 12, 13, 14)

# Groups for the X axis
coverage_groups = ((1, 3), (4, 6), (7, 9), (10, np.inf))

# Colors (copied from assigment specification)
plasma_colors = [(237/255, 252/255,  27/255),
                 (245/255, 135/255,  48/255),
                 (190/255,  48/255, 101/255),
                 ( 87/255,   0/255, 151/255),
                 ( 15/255,   0/255, 118/255)]

def make_colormap_single_channel(colors, channel, num_steps):
    """Calculate values for one RGB channel to make gradient colors."""
    num_colors = len(colors)
    # Number of gradient steps between each given color
    num_substeps = int((num_steps - num_colors) / (num_colors - 1))

    # Map from each color to next has num_substeps in between plus the edges
    submaps = [np.linspace(colors[i][channel], colors[i+1][channel],
                           num_substeps + 2)
               for i in range(num_colors - 1)]
    # Remove the duplicate internal boundary colors
    for i in range(num_colors - 2):
        submaps[i] = submaps[i][:-1]

    return np.concatenate(submaps, axis=None)

def make_colormap(colors, num_steps):
    """Calculate values to make gradient colors."""
    R = make_colormap_single_channel(colors, 0, num_steps)
    G = make_colormap_single_channel(colors, 1, num_steps)
    B = make_colormap_single_channel(colors, 2, num_steps)

    return np.vstack((R, G, B)).T

def dist_adj_factor(lims, panel_dim): return panel_dim / (lims[1] - lims[0])

def set_up_panel(figsize, loc, xlim, ylim):
    """Create a panel at (bottom, left, width, height)."""
    panel = plt.axes([loc[0] / figsize[0], loc[1] / figsize[1],
                     loc[2] / figsize[0], loc[3] / figsize[1]])
    panel.set_xlim(xlim)
    panel.set_ylim(ylim)
    return panel

def get_possible_x_pos(x_pos, width, min_dist, x_adj):
    """Get possible X positions for swarmplot."""
    # Step size for x moves
    increment = min_dist / 10 * x_adj
    # Order x shifts by distance from x_pos in ranage [-width/2, width/2]
    possible_positions = [x_pos]
    for shift in np.arange(increment, width / 2, increment):
        possible_positions.append(x_pos + shift)
        possible_positions.append(x_pos - shift)

    return possible_positions

def swarmplot(panel, points, x_pos, width, min_dist, x_adj, y_adj, gradient):
    """Make a swarmplot."""

    x_options = get_possible_x_pos(x_pos, width, min_dist, x_adj)
    # Tracking plotted points by Y bucket
    plotted_points, safe_x = dict(), dict()

    for point in points:
        y1, color = point
        # Buckets are half the minimum y-distance apart
        # Using the full minimum results in holes; skip too many x pos
        cur_bucket = int(y1 / min_dist * y_adj * 2)

        # Build up a list of positions sufficiently nearby to overlap
        nearby_points = []
        for nearby_buckets in range(cur_bucket - 2, cur_bucket + 3):
            if nearby_buckets not in plotted_points:
                plotted_points[nearby_buckets] = []
            nearby_points += plotted_points[nearby_buckets]
        
        plotted = False
        # Iterate starting at the first safe x position
        for i in range(safe_x.get(cur_bucket, 0), len(x_options)):
            x1 = x_options[i]

            # Check if using this x value would overlap with any nearby points
            has_overlap = False
            for x2, y2, in nearby_points:
                x_dist = (x2 - x1) * x_adj
                y_dist = (y2 - y1) * y_adj
                if (x_dist**2 + y_dist**2)**0.5 < min_dist:
                    has_overlap = True
                    # Short-circuit since this x position is impossible
                    break
            
            # Points are not allowed to overlap
            if not has_overlap:
                # Plot the point
                panel.plot(x1, y1, **POINT, color=gradient[color])
                plotted_points[cur_bucket].append((x1, y1))
                plotted = True

                # Other points in this bucket can't be plotted earlier
                safe_x[cur_bucket] = i
                # Short-curcuit since this point is plotted
                break
        
        # End early if this point failed to find a spot
        if not plotted:
            # Count how many points we managed to plot
            total_plotted = 0
            for bucketed_points in plotted_points.values():
                total_plotted += len(bucketed_points)
                
            print(f'{len(points) - total_plotted} points could not '
                    f'be plotted at position {x_pos}')
            break

def plot_main_panel(panel, swarm_point_max, reads, gradient, groups,
                    min_dist, x_adj, y_adj):
    """Plot left panel with cell types."""
    
    # Sort points into x-axis groups by subread coverage
    point_groups = {lims: [] for lims in groups}
    for read in reads.values():
        for lims in groups:
            if lims[0] <= read[1] <= lims[1]:
                point_groups[lims].append((read[0], read[2]))
                break
    
    # X ticks are built group-by-group
    ticks = {'pos': [], 'labels': []}
    for i in range(len(groups)):
        # Data coordinate for this x-group
        x = (i + 1) * 2 - 1
        # Only plot a maximum number of points
        points = random.sample(point_groups[groups[i]], swarm_point_max)
        swarmplot(panel, points, x, 1.8, min_dist, x_adj, y_adj, gradient)
        ticks['pos'].append(x)

        # Build a tick label for this range
        if groups[i][1] == np.inf:
            label = f'>={groups[i][0]}'
        else:
            label = f'{groups[i][0]}-{groups[i][1]}'
        ticks['labels'].append(label)

    panel.set_xticks(ticks['pos'], ticks['labels'])
    panel.set_xlabel('Subread Coverage')
    panel.set_ylabel('Identity (%)')

## == Main == ##

if __name__ == '__main__':
    identity, coverage, quality, output = get_params_command_line()

    reads = get_data(identity, coverage, quality, main_panel_lim['y'][0],
                     right_panel_ticks[0], right_panel_ticks[-1], NUM_STEPS)
    plasma = make_colormap(plasma_colors, NUM_STEPS)
    # Different axis limits/ranges mean the adjustment factor is different
    x_adj = dist_adj_factor(main_panel_lim['x'], main_panel_loc[2])
    y_adj = dist_adj_factor(main_panel_lim['y'], main_panel_loc[3])

    fig = plt.figure(figsize=figure_size, dpi=figure_dpi)
    main_panel = set_up_panel(figure_size, main_panel_loc,
                              main_panel_lim['x'], main_panel_lim['y'])
    right_panel = set_up_panel(figure_size, right_panel_loc,
                               right_panel_lim, right_panel_lim)

    plot_main_panel(main_panel, SWARM_SIZE, reads, plasma, coverage_groups,
                    POINT_OVERLAP_DIST, x_adj, y_adj)
    plot_right_panel(right_panel, right_panel_lim, NUM_STEPS, plasma, 
                     right_panel_ticks)

    fig.savefig(output)