function val = function_u(x, y)
    % FUNCTION_U_EXACT 定义解析解 (用于误差分析)
    
    val = sin(pi*x) .* sin(pi*y);
end