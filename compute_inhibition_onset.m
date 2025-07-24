function onset_time = compute_inhibition_onset(spike_times, pre_window, post_window, bin_size, threshold)
    % Compute baseline firing rate
    baseline_spikes = spike_times(spike_times >= pre_window(1) & spike_times < pre_window(2));
    baseline_rate = numel(baseline_spikes) / (pre_window(2) - pre_window(1));

    % Time bins for post-stimulus
    edges = post_window(1):bin_size:post_window(2);
    counts = histcounts(spike_times, edges);
    rates = counts / bin_size;
    times = edges(1:end-1);  % start time of each bin

    % Find first bin where rate drops below threshold of baseline
    drop_idx = find(rates < baseline_rate * threshold, 1);

    if ~isempty(drop_idx)
        onset_time = times(drop_idx);
    else
        onset_time = NaN;
    end
end
