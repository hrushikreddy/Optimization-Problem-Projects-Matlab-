function read_plot(data,depot)
    % uncomment next line for question 2.4:
    %data = data(:,1:2); % ONLY FIRST 2 COLS

    dims = size(data);
    m = dims(1);
    n = dims(2);
    
    xMat = data(:,1:2:end);
    yMat = data(:,2:2:end);

    title('Optimal Path')
    hold on
    plot(depot(1),depot(2),"*","color","blue")
    % SCATTER PLOT
    for i = 1:m
        scatter(xMat(i,:),yMat(i,:),"bla",".") % plot each row of the data as a scatter plot
        b = num2str(i); c = cellstr(b); % change number i to string
        dx = 0.5; dy = 0.5; % displacement so the text does not overlap with the data points
        text(xMat(i,:)+dx, yMat(i,:)+dy, c); % Assign a different number for each item
    end
    hold off

end