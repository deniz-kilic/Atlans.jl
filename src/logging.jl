# Taken from: https://github.com/Deltares/Wflow.jl
"""
    parse_loglevel(input_level::AbstractString)::LogLevel
    parse_loglevel(input_level::Integer)::LogLevel
Parse a log level from either an integer or string.
# Examples
    parse_loglevel("info") -> Logging.Info
    parse_loglevel(0) -> LogLevel(0) (== Logging.Info)
"""
function parse_loglevel(input_level::AbstractString)::LogLevel
    level = lowercase(input_level)
    if level == "debug"
        return Logging.Debug
    elseif level == "info"
        return Logging.Info
    elseif level == "warn"
        return Logging.Warn
    elseif level == "error"
        return Logging.Error
    else
        error("loglevel $input_level not recognized")
    end
end

parse_loglevel(input_level::Integer) = LogLevel(input_level)

function init_logger(loglevel, path_log, silent)
    loglevel = parse_loglevel(loglevel)
    log_handle = open(path_log, "w")
    # use ConsoleLogger instead of FileLogger to be able to print full stacktraces there
    # https://github.com/JuliaLogging/LoggingExtras.jl/issues/46#issuecomment-803716480
    ConsoleLogger(
        IOContext(log_handle, :compact => false, :limit => false, :color => false),
        loglevel,
        show_limited = false,
    )

    terminal_logger = if silent
        NullLogger()
    else
        # avoid debug information overflowing the terminal, -1 is the level of ProgressLogging
        # avoid double stacktraces by filtering the @error catch in Wflow.run
        EarlyFilteredLogger(
            log -> log.id !== :wflow_run,
            MinLevelLogger(TerminalLogger(), LogLevel(-1)),
        )
    end

    logger = TeeLogger(
        # avoid progresslogger and NCDatasets debug messages in the file
        EarlyFilteredLogger(
            log -> log.group !== :ProgressLogging && log._module != NCDatasets,
            terminal_logger,
        ),
        terminal_logger,
    )
    with_logger(logger) do
        @debug "Logger initialized." loglevel silent path_log
    end
    return logger, log_handle
end
