#function to calculate MA via algorithm based on a command
function run_with_timeout(command::Cmd, timeout::Float64)

    stdout_buffer = IOBuffer()
    stderr_buffer = IOBuffer()

    p = run(pipeline(command; stdout = stdout_buffer, stderr = stderr_buffer), wait = false)

    t_start = time()
    # check if output is returned or time limit is exceeded
    while isempty(readchomp(seekstart(stdout_buffer))) && (time() - t_start < timeout)
        sleep(0.1)
    end

    if isopen(p)
        kill(p, Base.SIGINT)
    end

	println(stdout_buffer.data)
	println(stderr_buffer.data)

	stdout_str = deepcopy(String(stdout_buffer.data))
	stderr_str = deepcopy(String(stderr_buffer.data))

    return stdout_str, stderr_str
end

#test the function
@show cmd = `./main --verbose nadph.mol` #34
@show stdout, stderr = run_with_timeout(cmd, 10.0)
println("stdout: " * stdout)
println("stderr: " * stderr)

# @info "Calculating MA values using molecular assembly algorithm..."
# @showprogress for i in nrow(complexities):-1:1
#     if ismissing(complexities.ma[i])
#         println(complexities.id[i])
#         cmd = `./main molfiles/$(complexities.id[i]).mol"`
#         stdout, stderr = run_with_timeout(cmd, 10.0)
#     end
# end
# @info "Done calculating MA values"
#
# open("data/complexities_$pathway.json", "w") do io
#     JSON3.write(io, complexities)
# end
