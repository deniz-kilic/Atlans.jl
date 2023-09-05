function nanpercentile(array, p)
    percentile(
        filter(!isnan, vec(array)),
        p
    )
end