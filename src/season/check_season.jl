# # 符合条件的两个相邻(t_diff <= 5) seasons, 合并endtime
# function LeftCombine_season(y_peak, y_end, len, date_beg, date_peak, date_end, 
#     i::Integer)
#     # ! this function is questionable
#     y_end[i + 1] = y_end[i];
#     date_end[i + 1] = date_end[i];
#     len[i + 1] = date_end[i + 1] - date_beg[i + 1] + 1;
#     if (y_peak[i] > y_peak[i + 1])
#         date_peak[i + 1] = date_peak[i];
#         y_peak[i + 1] = y_peak[i];
#     end
#     y_peak[i] = -9999.0; # flag
# end


function seasons_union(y_peak, y_end, len, date_beg, date_peak, date_end, 
    i::Integer)
    
    y_end[i] = y_end[i + 1];
    date_end[i] = date_end[i + 1];
    len[i] = date_end[i] - date_beg[i] + 1;

    if (y_peak[i] < y_peak[i + 1])
        date_peak[i] = date_peak[i + 1];
        y_peak[i] = y_peak[i + 1];
    end
    y_peak[i + 1] = -9999.0; ## flag
end


function check_season!(d; 
    rm_closed = true, rtrough_max = 0.7, r_min = 0.02)

    date_beg  = d[!, "beg"];
    date_end  = d[!, "end"];
    date_peak = d[!, "peak"];
    len       = d[!, "len"];
    y_beg     = d[!, "y_beg"];
    y_end     = d[!, "y_end"];
    y_peak    = d[!, "y_peak"];

    n = d.nrow();
    newdate;
    newval;
    
    for i = 1:n 
        t_diff = date_beg[i + 1] - date_end[i]; # 滞后时间，生长周期空窗期
        T2s    = date_end[i + 1] - date_beg[i]; # 两个生长周期长度

        T1_minVal = min(y_beg[i], y_end[i]);
        T2_minVal = min(y_beg[i + 1], y_end[i + 1]);
        maxY      = max(y_peak[i], y_peak[i + 1]);
        minY      = min(T1_minVal, T2_minVal);
        A         = maxY - minY;

        T1_h_right = y_peak[i] - y_end[i];         # T1, height of the right hand
        T2_h_left  = y_peak[i + 1] - y_beg[i + 1]; # T2, height of left hand

        trs  = A * r_min;
        trs2 = A * (r_min + 0.1);
        # # 1. y_end[i] >= trs || y_beg[i+1] >= trs 定义为高trough
        # #  diff =  abs(y_end[i] - y_beg[i+1]);
        # # bool is_HighTrough = y_end[i] >= trs || y_beg[i+1] >= trs;
        is_PreEndSmaller = (y_end[i] <= y_beg[i + 1]); # previous y_end smaller ?
        
        ## 现已考虑的情景包含如下：
        # 1. growing season日期交差                                  : t_diff < 0
        # 2. 两相邻val_troughs挨得太近（且中间不存在较大的peaks）    : t_diff <= 50
        # 3(temp). 两相邻val_troughs值相差过大; 且不重叠，则进行融合
        
        if (t_diff < 0) ##|| (is_HighTrough && delta_days < 0)
            # 如果两日期出现交叉，则重新洗牌`y_end[i]` and `y_beg[i+1]`
            # TODO: 补充出现这种情况的情景
            newdate = is_PreEndSmaller ? date_end[i] : date_beg[i + 1];
            newval = is_PreEndSmaller ? y_end[i] : y_beg[i + 1];

            date_end[i] = newdate; 
            date_beg[i + 1] = newdate;
            y_end[i] = newval;
            y_beg[i + 1] = newval;
            # # len[i] = date_end[i] - date_beg[i] + 1;
            # # len[i + 1] = date_end[i + 1] - date_beg[i + 1] + 1;
        end

        ## 坡脚修复大法的功能已经在`seasons_union`中得到了体现
        # IT-BCi, ending of 2005; 必须是相邻troughs
        if (t_diff > 0 && t_diff <= 150 &&
            T2_h_left <= trs && y_end[i] < y_beg[i + 1]) # 
            if (verbose); println("date_beg: from $(date_beg[i]) to $(date_beg[i+1])"); end
            # ? 为何添加 `y_end[i] < y_beg[i + 1]`的限制条件？
            # 向左拱进
            y_beg[i + 1] = y_beg[i];
            date_beg[i + 1] = date_beg[i];
            # # len[i + 1] = date_end[i + 1] - date_beg[i + 1] + 1;
        end

        # 处理相邻的peaks，进行融合操作
        # 1. 第二周期，左侧trough过高，超过`0.7*A` (BUG found)
        # 2. 第二周期，左侧r_min过小, （其ypeak > r_min + 0.1）
        # 3. 第一周期，右侧r_min过小,
        # 3. (not used) 如果相邻生长季返青日期过短，则执行合并操作
        # 这里用max，而非min，意在保护A较小的生长季
        if (rm_closed)
            is_closed = t_diff >= 0 && t_diff <= 150; ## 相差在5天内认为相邻
            
            # diff_right = y_peak[i + 1] - y_end[i + 1];
            # Rprintf("T1_h_right = %#4.2f, T2_h_left = % 4.2f, trs = %4.2f\n", T1_h_right, T2_h_left, trs);
            # a. 向左移动：如果T1_h_right过小，T1_h_left够大
            # b. 向左移动：如果T2_h_left过小，T2_h_right够大
            con_left = (T1_h_right <= trs && y_beg[i] < y_beg[i + 1] && (y_peak[i] - T1_minVal > trs2));
            con_right = (T2_h_left <= trs && y_end[i] > y_end[i + 1] && (y_peak[i + 1] - T1_minVal > trs2));
            
            if (is_closed && T2s <= 650 && 
                ((y_end[i] >= A * rtrough_max + T1_minVal) || con_right || con_left))
                # Rcout << i+1 << "正在进行融合" <<  endl;
                seasons_union(y_peak, y_end, len, date_beg, date_peak, date_end, i);
                # i++; ## 如果进行了融合，则跳过下一生长周期
                i += 1
                continue;
            end
            # if (is_closed && con_left)
            # {
            #     LeftCombine_season(y_peak, y_end, len, date_beg, date_peak, date_end, i);
            #     i++; # 如果进行了融合，则跳过下一生长周期
            #     continue;
            # }
        end
    end
    d[!, "len"] = date_end - date_beg + 1;
    # # CharacterVector names = d.attr("names");
    # # return max_dte;
end

export check_season!
