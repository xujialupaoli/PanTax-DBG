import sys

in_file = sys.argv[1]

with open(in_file,"r") as f:
    f.readline()
    for line in f:
        tokens = line.strip().split(":")
        if tokens[0].startswith("User"):
            cpu_user = float(tokens[-1])
        elif tokens[0].startswith("System"):
            cpu_system = float(tokens[-1])
        elif tokens[0].startswith("Elapsed"):
            elapsed_time = line.strip().replace("Elapsed (wall clock) time (h:mm:ss or m:ss): ", "")
            elapsed_time_parts = elapsed_time.strip().split(":")
            elapsed_in_seconds = float(elapsed_time_parts[-1])
            if len(elapsed_time_parts) > 1:
                elapsed_in_seconds += float(elapsed_time_parts[-2]) * 60
            if len(elapsed_time_parts) > 2:
                elapsed_in_seconds += float(elapsed_time_parts[-3]) * 3600

        elif tokens[0].startswith("Maximum"):
            memory = round(float(tokens[-1])/ 1024 / 1024, 1)
cpu_time = round(cpu_user + cpu_system, 1)

print('cpu(s) elapsed(s) memory(GB)\n')
print(str(cpu_time) + ' & ' + str('%.1f'% elapsed_in_seconds) + ' & ' + str(memory))

cpu_time = round(cpu_time/60/60, 1)
elapsed_in_hours = round(elapsed_in_seconds/60/60, 1)
print('\n' + '#' * 50 + '\n')
print('cpu(h) elapsed(h) memory(GB)\n')
print(str(cpu_time) + ' & ' + str(elapsed_in_hours) + ' & ' + str(memory))
print('\n\n')