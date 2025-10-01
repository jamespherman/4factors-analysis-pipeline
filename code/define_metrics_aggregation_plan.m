%% define_metrics_aggregation_plan.m
%
%   Defines a plan for aggregating per-neuron metrics and PSTH data.
%   This function serves as the single source of truth for the
%   `aggregate_neuron_metrics.m` script, specifying what data to extract
%   from the per-session files and how to structure it.
%
function metrics_plan = define_metrics_aggregation_plan()

    % Part 1: Define scalar metrics to be collected for every neuron.
    % Each entry defines a column in the output metrics table.
    metrics_plan.per_neuron_metrics(1).ColumnName = 'BaselineFR';
    metrics_plan.per_neuron_metrics(1).SourcePath = 'analysis.metrics.baseline_frs';
    
    metrics_plan.per_neuron_metrics(2).ColumnName = 'WaveformDuration';
    metrics_plan.per_neuron_metrics(2).SourcePath = ...
        'analysis.metrics.wf_metrics.peak_trough_ms';
    
    metrics_plan.per_neuron_metrics(3).ColumnName = 'IsSelected';
    metrics_plan.per_neuron_metrics(3).SourcePath = ...
        'analysis.selected_neurons';

    % Part 2: Define plan for aggregating time-resolved data (PSTHs) for
    % only the "selected" neurons.
    metrics_plan.psth_aggregation.SourcePath = ...
        'analysis.core_data.spikes';
    metrics_plan.psth_aggregation.Events = {'fixOn', 'targetOn', ...
        'fixOff', 'saccadeOnset', 'reward'};
    metrics_plan.psth_aggregation.DataField = 'rates';
    metrics_plan.psth_aggregation.SelectorPath = ...
        'analysis.selected_neurons';

end