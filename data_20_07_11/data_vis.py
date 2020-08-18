import csv
import matplotlib.pyplot as plt
import datetime

def read_data(file_name:str):
    with open(file_name, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter=',')
        data = list(reader)
    return data

def cumulative_incidences_per_province(data:list, province_names:list, days_with_zero_cases=0, plot=False):
    time = [i for i in range(days_with_zero_cases + 1)]
    cumulative_data = (days_with_zero_cases + 1) * [0]
    for line in data:
        if line[1] in province_names:
            time.append(time[-1] + 1)
            cumulative_data.append(cumulative_data[-1] + int(line[2]))
    max_data = max(cumulative_data)
    cumulative_data = [d/max_data for d in cumulative_data]
    if plot:
        title = ' - '
        plt.figure()
        plt.plot(time, cumulative_data)
        plt.title(title.join(province_names))
        plt.show()
    return time, cumulative_data


def incidences_per_province(data:list, province_names:list, days_with_zero_cases=0, plot=False):
    time = [i for i in range(days_with_zero_cases + 1)]
    new_data = (days_with_zero_cases + 1) * [0]
    for line in data:
        if line[1] in province_names:
            time.append(time[-1] + 1)
            new_data.append(int(line[2]))
    max_data = max(new_data)
    new_data = [d/max_data for d in new_data]
    if plot:
        title = ' - '
        plt.figure()
        plt.plot(time, new_data)
        plt.title(title.join(province_names))
        plt.show()
    return time, new_data


def save_data(time, data, file_name):
    with open(file_name, 'w', newline='') as f:
        writer = csv.writer(f)
        for t, d in zip(time, data):
            writer.writerow([t, d])
    print('Saving succeeded..')



if __name__ == '__main__':
    file_name = 'Belgium Covid19 interactive dashboard_Cases_Tijdreeks_cleaned.csv'
    data = read_data(file_name)
    data_date = datetime.date(2020, 3, 1)
    start_date = datetime.date(2020, 1, 24) # https://ec.europa.eu/info/live-work-travel-eu/health/coronavirus-response/timeline-eu-action_en
    date_difference = data_date - start_date
    print(date_difference.days)

    
    # Vlaanderen
    t_fl, cum_data_fl = cumulative_incidences_per_province(data, ['Antwerpen', 'VlaamsBrabant', 'OostVlaanderen', 'Limburg', 'WestVlaanderen'], days_with_zero_cases=date_difference.days, plot=True)
    # t_fl, data_fl = incidences_per_province(data, ['Antwerpen', 'VlaamsBrabant', 'OostVlaanderen', 'Limburg', 'WestVlaanderen'], days_with_zero_cases=date_difference.days, plot=True)
    save_data(t_fl, cum_data_fl, 'cum_cases_flanders.csv')

    # Antwerpen
    t_ant, cum_data_ant = cumulative_incidences_per_province(data, ['Antwerpen'], days_with_zero_cases=date_difference.days, plot=True)
    save_data(t_ant, cum_data_ant, 'cum_cases_antwerp.csv')

    # Vlaams-Brabant
    t_vb, cum_data_vb = cumulative_incidences_per_province(data, ['VlaamsBrabant'], days_with_zero_cases=date_difference.days, plot=True)
    save_data(t_vb, cum_data_vb, 'cum_cases_flemBrab.csv')

    # Oost-Vlaanderen
    t_ef, cum_data_ef = cumulative_incidences_per_province(data, ['OostVlaanderen'], days_with_zero_cases=date_difference.days, plot=True)
    save_data(t_ef, cum_data_ef, 'cum_cases_eastFl.csv')

    # Limburg
    t_lim, cum_data_lim = cumulative_incidences_per_province(data, ['Limburg'], days_with_zero_cases=date_difference.days, plot=True)
    save_data(t_lim, cum_data_lim, 'cum_cases_limburg.csv')

    # West-Vlaanderen
    t_wv, cum_data_wv = cumulative_incidences_per_province(data, ['WestVlaanderen'], days_with_zero_cases=date_difference.days, plot=True)
    save_data(t_wv, cum_data_wv, 'cum_cases_westFl.csv')
