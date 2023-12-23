merge_segments<-function(segments) {
  merged_segments <- segments[1, ]  # Initialize merged segments with the first segment
  for (i in 2:(nrow(segments)-1)) { 
	curr_segment <- segments[i, ]
	next_segment <- segments[i+1, ]
    prev_segment <- merged_segments[nrow(merged_segments), ]
    if ((curr_segment$seg.mean > 1.5 && prev_segment$seg.mean > 1.5) || (curr_segment$seg.mean < 1.4 && prev_segment$seg.mean < 1.4) || (curr_segment$seg.mean < 1.5 && prev_segment$seg.mean > 1.5 && curr_segment$Length<5000 && next_segment$seg.mean>1.5) || (curr_segment$seg.mean > 1.5 && prev_segment$seg.mean < 1.5 && curr_segment$Length<5000 && next_segment$seg.mean<1.5)) {
      # Merge segments if both have seg.mean > 1.5, or separated by a num.mark<10 segment
      merged_segments[nrow(merged_segments), "loc.end"] <- curr_segment$loc.end
      merged_segments[nrow(merged_segments), "Length"] <- prev_segment$Length + curr_segment$Length
      merged_segments[nrow(merged_segments), "num.mark"] <- prev_segment$num.mark + curr_segment$num.mark
      merged_segments[nrow(merged_segments), "seg.mean"] <- (prev_segment$seg.mean * prev_segment$Length + curr_segment$seg.mean * curr_segment$Length) / (prev_segment$Length + curr_segment$Length)
    } else {
      # Append current segment as a new merged segment
      merged_segments <- rbind(merged_segments, curr_segment)
	} 
  }
i<- nrow(segments)
curr_segment <- merged_segments[nrow(merged_segments), ]
next_segment <- segments[i, ]
prev_segment <- merged_segments[nrow(merged_segments)-1, ]
if ((curr_segment$seg.mean > 1.5 && prev_segment$seg.mean > 1.5) || (curr_segment$seg.mean < 1.4 && prev_segment$seg.mean < 1.4) || (curr_segment$seg.mean < 1.5 && prev_segment$seg.mean > 1.5 && curr_segment$Length<5000 && next_segment$seg.mean>1.5) || (curr_segment$seg.mean > 1.5 && prev_segment$seg.mean < 1.5 && curr_segment$Length<5000 && next_segment$seg.mean<1.5)) {
  # Merge segments if both have seg.mean > 1.5, or separated by a num.mark<10 segment
  merged_segments[nrow(merged_segments), "loc.end"] <- curr_segment$loc.end
  merged_segments[nrow(merged_segments), "Length"] <- prev_segment$Length + curr_segment$Length
  merged_segments[nrow(merged_segments), "num.mark"] <- prev_segment$num.mark + curr_segment$num.mark
  merged_segments[nrow(merged_segments), "seg.mean"] <- (prev_segment$seg.mean * prev_segment$Length + curr_segment$seg.mean * curr_segment$Length) / (prev_segment$Length + curr_segment$Length)
  } else {
    # Append current segment as a new merged segment
    merged_segments <- rbind(merged_segments, next_segment)
  }
	return(merged_segments)
}
